"""FastAPI 入口：批次提交、批次结果查询、批次列表、批次风险对比。"""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import RedirectResponse

from . import config, db
from .schemas import (
    AugmentationPlanOut,
    AugmentationPlanRequest,
    BatchCompareResult,
    BatchCreate,
    BatchList,
    BatchResult,
    BatchSummary,
    PairOut,
    PoolPlanOut,
    PoolPlanRequest,
    PrimerOut,
    ReplacementRankingOut,
    ReplacementRankingRequest,
    ThresholdMapOut,
    ThresholdMapRequest,
)
from .services import (
    build_batch_threshold_map,
    compare_batches,
    plan_batch_augmentation,
    plan_batch_pools,
    rank_replacements,
    run_batch,
)


@asynccontextmanager
async def lifespan(_app: FastAPI):
    # 依赖与数据库初始化纳入容器启动流程
    db.init_db(reset=config.RESET_ON_START)
    yield


app = FastAPI(
    title="PCR 引物二聚体筛查 API",
    description=(
        "提交一批带名称的 IUPAC DNA 引物，枚举自二聚体与两两异二聚体的全部"
        "错位对齐，按 IUPAC 碱基集合判断互补，返回每对最危险对齐、连续互补长度、"
        "涉及 3′端的连续长度与风险原因。"
    ),
    version="1.0.0",
    lifespan=lifespan,
)


@app.get("/", include_in_schema=False)
def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")


@app.get("/health", tags=["meta"], summary="健康检查")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/api/batches",
    status_code=201,
    tags=["batches"],
    summary="提交引物批次并同步完成筛查",
    response_model=BatchResult,
)
def create_batch(payload: BatchCreate) -> BatchResult:
    if len(payload.primers) > config.MAX_PRIMERS_PER_BATCH:
        raise HTTPException(
            status_code=422,
            detail=(
                f"引物数量 {len(payload.primers)} 超过上限 "
                f"{config.MAX_PRIMERS_PER_BATCH}"
            ),
        )

    seen: set[str] = set()
    duplicates: list[str] = []
    primers: list[tuple[str, str]] = []
    for p in payload.primers:
        key = p.name
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        primers.append((p.name, p.sequence))
    if duplicates:
        raise HTTPException(
            status_code=422,
            detail=f"批次内引物名称重复: {', '.join(duplicates)}",
        )

    batch_id = run_batch(
        primers=primers,
        continuous_threshold=payload.continuous_threshold,
        three_prime_threshold=payload.three_prime_threshold,
    )
    return _load_result(batch_id, risky_only=False)  # type: ignore[return-value]


@app.get(
    "/api/batches",
    tags=["batches"],
    summary="列出历史批次（按创建时间倒序）",
    response_model=BatchList,
)
def list_batches(
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
) -> BatchList:
    with db.get_conn() as conn:
        total = db.count_batches(conn)
        rows = db.list_batches(conn, limit, offset)
    return BatchList(
        total=total,
        limit=limit,
        offset=offset,
        batches=[_summary(r) for r in rows],
    )


@app.get(
    "/api/batches/{batch_id}",
    tags=["batches"],
    summary="查询批次完整筛查结果",
    response_model=BatchResult,
    responses={404: {"description": "批次不存在"}},
)
def get_batch(
    batch_id: str,
    risky_only: bool = Query(
        False, description="为 true 时只返回 high/medium 风险配对"
    ),
) -> BatchResult:
    return _load_result(batch_id, risky_only=risky_only)


def _summary(row) -> BatchSummary:
    return BatchSummary(**dict(row))


@app.get(
    "/api/batches/{base_id}/compare/{target_id}",
    tags=["batches"],
    summary="对比两个批次的配对风险变化（base → target）",
    response_model=BatchCompareResult,
    responses={404: {"description": "任一批次不存在"}},
)
def compare_batch(base_id: str, target_id: str) -> BatchCompareResult:
    result = compare_batches(base_batch_id=base_id, target_batch_id=target_id)
    if result is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    return BatchCompareResult(**result)


@app.post(
    "/api/batches/{batch_id}/pool-plan",
    tags=["batches"],
    summary="基于已保存批次的异二聚体结果做引物分池",
    response_model=PoolPlanOut,
    responses={
        404: {"description": "批次不存在"},
        422: {"description": "参数非法或池容量总和不足"},
    },
)
def create_pool_plan(batch_id: str, payload: PoolPlanRequest) -> PoolPlanOut:
    try:
        result = plan_batch_pools(
            batch_id=batch_id,
            pool_count=payload.pool_count,
            pool_capacity=payload.pool_capacity,
            isolated_risk_levels=payload.isolated_risk_levels,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    if result is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    return PoolPlanOut(**result)


@app.post(
    "/api/batches/{batch_id}/replacement-ranking",
    tags=["batches"],
    summary="候选替换引物的内存模拟与非支配分层（不落库）",
    response_model=ReplacementRankingOut,
    responses={
        404: {"description": "批次不存在"},
        422: {"description": "目标引物不存在或候选标签重复"},
    },
)
def create_replacement_ranking(
    batch_id: str, payload: ReplacementRankingRequest
) -> ReplacementRankingOut:
    seen: set[str] = set()
    duplicates: list[str] = []
    candidates: list[tuple[str, str]] = []
    for c in payload.candidates:
        if c.label in seen:
            duplicates.append(c.label)
        seen.add(c.label)
        candidates.append((c.label, c.sequence))
    if duplicates:
        raise HTTPException(
            status_code=422,
            detail=f"候选标签重复: {', '.join(duplicates)}",
        )

    try:
        result = rank_replacements(
            batch_id=batch_id,
            target_name=payload.target_name,
            candidates=candidates,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    if result is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    return ReplacementRankingOut(**result)


@app.post(
    "/api/batches/{batch_id}/augmentation-plan",
    tags=["batches"],
    summary="候选引物增量筛选：冲突图上的确定性分支限界精确解（不落库）",
    response_model=AugmentationPlanOut,
    responses={
        404: {"description": "批次不存在"},
        422: {"description": "候选名称非法/重复或数量超限"},
    },
)
def create_augmentation_plan(
    batch_id: str, payload: AugmentationPlanRequest
) -> AugmentationPlanOut:
    seen: set[str] = set()
    duplicates: list[str] = []
    candidates: list[tuple[str, str]] = []
    for c in payload.candidates:
        if c.name in seen:
            duplicates.append(c.name)
        seen.add(c.name)
        candidates.append((c.name, c.sequence))
    if duplicates:
        raise HTTPException(
            status_code=422,
            detail=f"候选引物名称重复: {', '.join(duplicates)}",
        )

    try:
        result = plan_batch_augmentation(
            batch_id=batch_id,
            max_selected=payload.max_selected,
            candidates=candidates,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    if result is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    return AugmentationPlanOut(**result)


@app.post(
    "/api/batches/{batch_id}/threshold-map",
    tags=["batches"],
    summary="已保存批次的阈值扫描地图（二维差分聚合，不重算对齐）",
    response_model=ThresholdMapOut,
    responses={
        404: {"description": "批次不存在"},
        422: {"description": "阈值范围非法或组合数超限"},
    },
)
def create_threshold_map(
    batch_id: str, payload: ThresholdMapRequest
) -> ThresholdMapOut:
    result = build_batch_threshold_map(
        batch_id=batch_id,
        continuous_start=payload.continuous_start,
        continuous_end=payload.continuous_end,
        three_prime_start=payload.three_prime_start,
        three_prime_end=payload.three_prime_end,
    )
    if result is None:
        raise HTTPException(status_code=404, detail="批次不存在")
    return ThresholdMapOut(**result)


def _load_result(batch_id: str, *, risky_only: bool) -> BatchResult:
    with db.get_conn() as conn:
        batch = db.get_batch(conn, batch_id)
        if batch is None:
            raise HTTPException(status_code=404, detail="批次不存在")
        primers = db.get_primers(conn, batch_id)
        pairs = db.get_pairs(conn, batch_id, risky_only=risky_only)

    return BatchResult(
        **dict(batch),
        primers=[PrimerOut(**dict(p)) for p in primers],
        pairs=[PairOut(**db.serialize_pair(r)) for r in pairs],
    )
