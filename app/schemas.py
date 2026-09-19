"""Pydantic 请求/响应模型与输入校验。"""

from __future__ import annotations

from pydantic import BaseModel, Field, field_validator, model_validator

from . import config
from .iupac import validate_iupac


class PrimerIn(BaseModel):
    name: str = Field(..., description="引物名称（非空，批次内唯一）")
    sequence: str = Field(
        ..., description="5'→3' IUPAC DNA 序列（允许 ATCGRYSWKMBDHVN）"
    )

    @field_validator("name")
    @classmethod
    def _name(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("引物名称不能为空")
        if len(v) > 100:
            raise ValueError("引物名称最长 100 字符")
        return v

    @field_validator("sequence")
    @classmethod
    def _sequence(cls, v: str) -> str:
        seq = validate_iupac(v)
        if len(seq) > config.MAX_PRIMER_LENGTH:
            raise ValueError(
                f"引物长度 {len(seq)} 超过上限 {config.MAX_PRIMER_LENGTH}"
            )
        return seq


class BatchCreate(BaseModel):
    primers: list[PrimerIn] = Field(
        ..., min_length=1, description="带名称的 IUPAC DNA 引物批次"
    )
    continuous_threshold: int = Field(
        ...,
        ge=1,
        le=config.MAX_PRIMER_LENGTH,
        description="连续互补长度阈值（达到即提示结合稳定风险）",
    )
    three_prime_threshold: int = Field(
        ...,
        ge=1,
        le=config.MAX_PRIMER_LENGTH,
        description="3′端连续互补长度阈值（达到即提示延伸风险，风险等级 high）",
    )


class AlignmentOut(BaseModel):
    top: str
    match: str
    bottom: str


class PairOut(BaseModel):
    pair_index: int
    name_top: str
    name_bottom: str
    is_self_dimer: bool
    alignments_checked: int
    best_shift: int = Field(
        description=(
            "最危险对齐的带符号错位：bottom 3′端碱基的全局列坐标（top 5′端 "
            "T0 位于列 0）。0=齐平，>0=bottom 3′端伸入 top，<0=悬空于 top 5′左侧"
        )
    )
    overlap_len: int
    continuous_complementary_length: int
    three_prime_run_length: int
    top3_run_length: int
    bottom3_run_length: int
    involved_three_prime_ends: list[str]
    risk_level: str = Field(description="high / medium / low")
    risk_reasons: list[str]
    alignment: AlignmentOut


class PrimerOut(BaseModel):
    position: int
    name: str
    sequence: str


class BatchSummary(BaseModel):
    id: str
    created_at: str
    continuous_threshold: int
    three_prime_threshold: int
    primer_count: int
    pair_count: int
    high_risk_count: int
    medium_risk_count: int
    low_risk_count: int


class BatchResult(BatchSummary):
    primers: list[PrimerOut]
    pairs: list[PairOut]


class BatchList(BaseModel):
    total: int
    limit: int
    offset: int
    batches: list[BatchSummary]


_RISK_LEVELS = ("high", "medium", "low")


class PoolPlanRequest(BaseModel):
    pool_count: int = Field(
        ..., ge=1, le=10000, description="池数量（至少 1）"
    )
    pool_capacity: int = Field(
        ..., ge=1, le=10000, description="每池容量（每池最多容纳的引物数）"
    )
    isolated_risk_levels: list[str] = Field(
        ...,
        min_length=1,
        description=(
            "需隔离的异二聚体风险等级（仅 high/medium/low）；达到这些等级的"
            "异二聚体两两构成冲突图边，自二聚体忽略。重复值自动去重"
        ),
    )

    @field_validator("isolated_risk_levels")
    @classmethod
    def _isolated_risk_levels(cls, values: list[str]) -> list[str]:
        normalized: list[str] = []
        bad: list[str] = []
        for value in values:
            level = value.strip().lower()
            if level not in _RISK_LEVELS:
                bad.append(value)
            elif level not in normalized:
                normalized.append(level)
        if bad:
            raise ValueError(
                "风险等级仅支持 high/medium/low，非法值: "
                + ", ".join(repr(v) for v in bad)
            )
        return normalized


class PoolConflictEdgeOut(BaseModel):
    pool: int = Field(description="冲突所在池号（0 起）")
    position_a: int = Field(description="引物 A 在批次中的位置")
    position_b: int = Field(description="引物 B 在批次中的位置")
    name_a: str
    name_b: str
    risk_level: str = Field(description="high / medium / low")


class PoolOut(BaseModel):
    pool: int = Field(description="池号（0 起）")
    capacity: int
    members: list[str] = Field(description="池内引物名称（按批次位置升序）")
    member_positions: list[int] = Field(
        description="与 members 一一对应的批次位置"
    )
    conflict_count: int = Field(description="该池内部的同池冲突边数")


class PoolPlanOut(BaseModel):
    batch_id: str
    pool_count: int
    pool_capacity: int
    isolated_risk_levels: list[str]
    primer_count: int
    assigned_count: int = Field(description="实际分配引物数（容量允许时等于引物数）")
    total_conflicts: int = Field(description="全部池的同池冲突边总数")
    pools: list[PoolOut]
    conflicts: list[PoolConflictEdgeOut]


class ComparePairOut(BaseModel):
    pair_key: str = Field(
        description="稳定配对键：自二聚体为单名，异二聚体为按名称字典序的 A|B"
    )
    name_a: str = Field(description="配对键第一名称（异二聚体为字典序较小者）")
    name_b: str = Field(description="配对键第二名称（自二聚体与 name_a 相同）")
    is_self_dimer: bool
    sequence_changed: bool = Field(
        description="配对中任一同名引物在两批次间序列不同"
    )
    base_risk_level: str | None = Field(
        description="base 批次风险等级（introduced 时为 null）"
    )
    target_risk_level: str | None = Field(
        description="target 批次风险等级（resolved 时为 null）"
    )


class CompareSummary(BaseModel):
    base_batch_id: str
    target_batch_id: str
    base_pair_count: int = Field(description="base 批次配对总数")
    target_pair_count: int = Field(description="target 批次配对总数")
    common_pair_count: int = Field(description="两批次共有的配对键数")
    unchanged_pair_count: int = Field(description="共有配对中风险等级不变的数量")
    escalated_count: int
    deescalated_count: int
    introduced_count: int
    resolved_count: int


class BatchCompareResult(BaseModel):
    summary: CompareSummary
    escalated: list[ComparePairOut] = Field(
        description="共有配对中风险等级上升（low<medium<high）的明细，按配对名排序"
    )
    deescalated: list[ComparePairOut] = Field(
        description="共有配对中风险等级下降的明细，按配对名排序"
    )
    introduced: list[ComparePairOut] = Field(
        description="仅存在于 target 批次的配对明细，按配对名排序"
    )
    resolved: list[ComparePairOut] = Field(
        description="仅存在于 base 批次的配对明细，按配对名排序"
    )


class ErrorOut(BaseModel):
    detail: str


class ReplacementCandidateIn(BaseModel):
    label: str = Field(
        ..., description="候选标签（非空，本次请求内唯一），仅用于回传标识"
    )
    sequence: str = Field(
        ..., description="5'→3' IUPAC DNA 序列（允许 ATCGRYSWKMBDHVN）"
    )

    @field_validator("label")
    @classmethod
    def _label(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("候选标签不能为空")
        if len(v) > 100:
            raise ValueError("候选标签最长 100 字符")
        return v

    @field_validator("sequence")
    @classmethod
    def _sequence(cls, v: str) -> str:
        seq = validate_iupac(v)
        if len(seq) > config.MAX_PRIMER_LENGTH:
            raise ValueError(
                f"候选序列长度 {len(seq)} 超过上限 {config.MAX_PRIMER_LENGTH}"
            )
        return seq


class ReplacementRankingRequest(BaseModel):
    target_name: str = Field(
        ..., description="待替换引物名称（必须是批次内已存在的引物名）"
    )
    candidates: list[ReplacementCandidateIn] = Field(
        ..., min_length=2, max_length=50, description="2～50 条带唯一标签的候选序列"
    )

    @field_validator("target_name")
    @classmethod
    def _target_name(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("待替换引物名称不能为空")
        return v


class ReplacementPairOut(PairOut):
    """替换模拟中受影响的配对（候选自二聚体或候选与其余引物的异二聚体）。"""


class ReplacementCandidateOut(BaseModel):
    label: str
    sequence: str
    layer: int = Field(description="非支配层级（1 起，越小越优）")
    rank: int = Field(description="最终名次（1 起；先按层级，再按输入顺序）")
    high_pair_count: int = Field(
        description="按批次原阈值判定为 high 的受影响配对数（含候选自二聚体）"
    )
    max_three_prime_run: int = Field(
        description="受影响配对中的最大 3′连续互补长度"
    )
    max_continuous_run: int = Field(
        description="受影响配对中的最长连续互补长度"
    )
    high_pairs: list[ReplacementPairOut] = Field(
        description="全部 high 受影响配对（按配对索引升序）"
    )
    max_three_prime_pairs: list[ReplacementPairOut] = Field(
        description="3′连续互补达到 max_three_prime_run 的配对（按配对索引升序）"
    )
    max_continuous_pairs: list[ReplacementPairOut] = Field(
        description="连续互补达到 max_continuous_run 的配对（按配对索引升序）"
    )


class ReplacementRankingOut(BaseModel):
    batch_id: str
    target_name: str
    continuous_threshold: int = Field(description="批次原连续互补阈值")
    three_prime_threshold: int = Field(description="批次原 3′连续互补阈值")
    primer_count: int = Field(description="批次引物数（模拟中保持不变）")
    candidate_count: int
    candidates: list[ReplacementCandidateOut] = Field(
        description="候选排名（按层、层内按请求输入顺序稳定排序）"
    )


class AugmentationCandidateIn(PrimerIn):
    """增量候选：复用 PrimerIn 的名称与 IUPAC 序列校验。"""


class AugmentationPlanRequest(BaseModel):
    max_selected: int = Field(
        ..., ge=1, le=18, description="最多接纳（入选）的候选数（1～18）"
    )
    candidates: list[AugmentationCandidateIn] = Field(
        ..., min_length=2, max_length=18, description="2～18 条名称唯一的候选引物"
    )


class AugmentationBlockingPairOut(PairOut):
    """阻断候选入选的一条 high/medium 配对。"""

    counterpart_kind: str = Field(
        description="配对另一端的类别：self=候选自二聚体，batch=原批次引物，candidate=入选候选"
    )
    counterpart_name: str = Field(
        description="配对另一端名称（自二聚体时与候选同名）"
    )


class AugmentationSelectedOut(BaseModel):
    position: int = Field(description="候选在请求中的序号（0 起）")
    name: str
    sequence: str


class AugmentationRejectedOut(BaseModel):
    position: int = Field(description="候选在请求中的序号（0 起）")
    name: str
    sequence: str
    reason: str = Field(
        description="淘汰原因：risk=high/medium 配对阻断；capacity=无风险冲突但接纳数已满"
    )
    limited_by_capacity: bool = Field(
        description="true 表示候选本身无阻断配对，仅因不超过 max_selected 的数量上限未入选"
    )
    blocking_pairs: list[AugmentationBlockingPairOut] = Field(
        description=(
            "risk 淘汰时的全部阻断配对（候选自二聚体、与原批次引物、与入选候选），"
            "按全局配对顺序排列；capacity 淘汰时为空列表"
        )
    )


class AugmentationPlanOut(BaseModel):
    batch_id: str
    continuous_threshold: int = Field(description="沿用的批次原连续互补阈值")
    three_prime_threshold: int = Field(description="沿用的批次原 3′端连续互补阈值")
    max_selected: int = Field(description="请求给定的最多接纳数")
    candidate_count: int
    selected_count: int
    rejected_count: int
    selected: list[AugmentationSelectedOut] = Field(
        description="入选候选（按请求顺序）"
    )
    rejected: list[AugmentationRejectedOut] = Field(
        description="未入选候选（按请求顺序）"
    )


class ThresholdMapRequest(BaseModel):
    continuous_start: int = Field(
        ...,
        ge=1,
        le=config.THRESHOLD_MAP_MAX_THRESHOLD,
        description="连续互补阈值起点（含，1～200）",
    )
    continuous_end: int = Field(
        ...,
        ge=1,
        le=config.THRESHOLD_MAP_MAX_THRESHOLD,
        description="连续互补阈值终点（含，1～200）",
    )
    three_prime_start: int = Field(
        ...,
        ge=1,
        le=config.THRESHOLD_MAP_MAX_THRESHOLD,
        description="3′端连续互补阈值起点（含，1～200）",
    )
    three_prime_end: int = Field(
        ...,
        ge=1,
        le=config.THRESHOLD_MAP_MAX_THRESHOLD,
        description="3′端连续互补阈值终点（含，1～200）",
    )

    @model_validator(mode="after")
    def _check_ranges(self) -> "ThresholdMapRequest":
        if self.continuous_start > self.continuous_end:
            raise ValueError("continuous_start 不能大于 continuous_end")
        if self.three_prime_start > self.three_prime_end:
            raise ValueError("three_prime_start 不能大于 three_prime_end")
        cells = (self.continuous_end - self.continuous_start + 1) * (
            self.three_prime_end - self.three_prime_start + 1
        )
        if cells > config.THRESHOLD_MAP_MAX_CELLS:
            raise ValueError(
                f"阈值组合数 {cells} 超过上限 {config.THRESHOLD_MAP_MAX_CELLS}"
            )
        return self


class ThresholdMapCellOut(BaseModel):
    continuous_threshold: int = Field(description="该单元的连续互补阈值")
    three_prime_threshold: int = Field(description="该单元的 3′端连续互补阈值")
    high_count: int = Field(description="该阈值组合下判定为 high 的配对数")
    medium_count: int = Field(description="该阈值组合下判定为 medium 的配对数")
    low_count: int = Field(description="该阈值组合下判定为 low 的配对数")
    is_batch_threshold: bool = Field(
        description="该单元是否为批次提交时使用的原始阈值组合"
    )


class ThresholdMapOut(BaseModel):
    batch_id: str
    continuous_start: int
    continuous_end: int
    three_prime_start: int
    three_prime_end: int
    cell_count: int = Field(description="阈值组合（单元）总数")
    pair_count: int = Field(description="批次配对总数（每单元三级计数之和）")
    batch_continuous_threshold: int = Field(description="批次原连续互补阈值")
    batch_three_prime_threshold: int = Field(description="批次原 3′端连续互补阈值")
    cells: list[ThresholdMapCellOut] = Field(
        description="按 3′端阈值、连续互补阈值升序的逐单元三级计数"
    )
