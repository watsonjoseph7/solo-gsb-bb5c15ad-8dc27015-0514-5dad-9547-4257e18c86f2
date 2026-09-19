"""分池接口端到端测试与 DSatur 贪心算法单测。"""

from __future__ import annotations

import os

import pytest

from app.pooling import build_pool_plan


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("pool_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},   # pos 0
        {"name": "R", "sequence": "GGGGGTTTTT"},   # pos 1
        {"name": "X", "sequence": "AAAAAAAAAA"},   # pos 2
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}
# 异二聚体风险：F-R=high、R-X=high、F-X=low


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _plan(client, batch_id, **over):
    body = {
        "pool_count": 2,
        "pool_capacity": 2,
        "isolated_risk_levels": ["high"],
    }
    body.update(over)
    return client.post(f"/api/batches/{batch_id}/pool-plan", json=body)


# ---------- 正常分池 ----------

def test_separates_high_risk_pair_and_assigns_all(client, batch_id):
    data = _plan(client, batch_id).json()

    assert data["primer_count"] == 3
    assert data["assigned_count"] == 3
    assert data["total_conflicts"] == 0
    assert data["conflicts"] == []

    pool_of = {m: p["pool"] for p in data["pools"] for m in p["members"]}
    assert pool_of == {"F": 1, "R": 0, "X": 1}  # 首批 R 单独占池 0

    # 各池成员按批次位置升序
    for pool in data["pools"]:
        assert pool["member_positions"] == sorted(pool["member_positions"])
        assert len(pool["members"]) == len(pool["member_positions"])
        assert len(pool["members"]) <= pool["capacity"]

    # 容量总和（4）允许时必须全部分配
    assert sum(len(p["members"]) for p in data["pools"]) == 3


def test_low_risk_level_not_isolated_gives_no_conflicts(client, batch_id):
    # 仅隔离 medium：图中没有 medium 异二聚体，冲突图无边
    data = _plan(client, batch_id, isolated_risk_levels=["medium"]).json()
    assert data["total_conflicts"] == 0
    assert data["assigned_count"] == 3
    assert sum(p["conflict_count"] for p in data["pools"]) == 0


def test_low_edge_becomes_conflict_when_isolated(client, batch_id):
    # 隔离全部等级：F-X 的 low 边也进入冲突图（完整 K3）
    data = _plan(
        client, batch_id, isolated_risk_levels=["high", "medium", "low"]
    ).json()
    assert data["assigned_count"] == 3
    # 2 池容量 2 装 K3 必有且仅有一条同池边
    assert data["total_conflicts"] == 1
    edge = data["conflicts"][0]
    assert edge["name_a"] in {"F", "R", "X"}
    # 两端确实在 edge.pool 所指的同一池中
    edge_pool = next(p for p in data["pools"] if p["pool"] == edge["pool"])
    assert edge["name_a"] in edge_pool["members"]
    assert edge["name_b"] in edge_pool["members"]
    # 其余池不含这条边的任一端
    for pool in data["pools"]:
        if pool["pool"] != edge["pool"]:
            assert edge["name_a"] not in pool["members"]
            assert edge["name_b"] not in pool["members"]
    # 池冲突计数之和等于总冲突数
    assert sum(p["conflict_count"] for p in data["pools"]) == 1


# ---------- 自二聚体被忽略 ----------

def test_self_dimers_are_ignored(client):
    # AAAATTTT 的自二聚体为 high；单引物批次不应产生任何冲突
    bid = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "S", "sequence": "AAAATTTT"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 3,
        },
    ).json()["id"]
    data = _plan(client, bid, pool_count=1, pool_capacity=1).json()
    assert data["primer_count"] == 1
    assert data["assigned_count"] == 1
    assert data["total_conflicts"] == 0
    assert data["pools"][0]["members"] == ["S"]


# ---------- 错误码 ----------

def test_batch_not_found_404(client):
    resp = _plan(client, "no-such-batch")
    assert resp.status_code == 404
    assert "不存在" in resp.json()["detail"]


def test_capacity_shortage_422(client, batch_id):
    resp = _plan(client, batch_id, pool_count=1, pool_capacity=2)
    assert resp.status_code == 422
    assert "不足" in resp.json()["detail"]


def test_exact_capacity_ok(client, batch_id):
    # 容量恰好等于引物数：合法，必须全部分配
    data = _plan(client, batch_id, pool_count=3, pool_capacity=1).json()
    assert data["assigned_count"] == 3
    assert all(len(p["members"]) <= 1 for p in data["pools"])
    assert data["total_conflicts"] == 0


def test_invalid_risk_level_422(client, batch_id):
    resp = _plan(client, batch_id, isolated_risk_levels=["critical"])
    assert resp.status_code == 422


def test_invalid_pool_params_422(client, batch_id):
    assert _plan(client, batch_id, pool_count=0).status_code == 422
    assert _plan(client, batch_id, pool_capacity=0).status_code == 422
    assert _plan(client, batch_id, isolated_risk_levels=[]).status_code == 422


# ---------- 稳定性 ----------

def test_same_input_stable(client, batch_id):
    a = _plan(client, batch_id).json()
    b = _plan(client, batch_id).json()
    assert a["pools"] == b["pools"]
    assert a["conflicts"] == b["conflicts"]
    assert a["total_conflicts"] == b["total_conflicts"]


def test_risk_levels_dedup_preserves_order(client, batch_id):
    data = _plan(
        client, batch_id, isolated_risk_levels=["medium", "high", "medium"]
    ).json()
    assert data["isolated_risk_levels"] == ["medium", "high"]


# ---------- DSatur 贪心算法单测 ----------

def _assignments(plan):
    return [plan.assignments[i] for i in range(len(plan.assignments))]


def test_dsatur_triangle_forced_one_conflict():
    # 完整三角形，2 池容量 2：必有 1 条冲突；位置 0 最先入池 0
    plan = build_pool_plan(
        primer_count=3,
        edges=[(0, 1, "high"), (1, 2, "high"), (0, 2, "high")],
        pool_count=2,
        pool_capacity=2,
    )
    assert plan.total_conflicts == 1
    assert _assignments(plan)[0] == 0
    assert len(plan.conflicts) == 1
    pool_id, edge = plan.conflicts[0]
    assert edge.risk_level == "high"
    assert plan.assignments[edge.position_a] == plan.assignments[edge.position_b] == pool_id


def test_dsatur_saturation_drives_next_choice():
    # 三角形 {0,1,2} + 尾巴 1-3；4 个顶点，3 池容量 2，理论可 0 冲突。
    # 第二选必须是与 0 相邻且度数更大的 1（而非位置更靠前的 0 邻居 2），
    # 且已见邻池数主导：最终 K3 三色分开、尾巴随 0 同池。
    edges = [(0, 1, "high"), (1, 2, "high"), (0, 2, "high"), (1, 3, "medium")]
    plan = build_pool_plan(
        primer_count=4, edges=edges, pool_count=3, pool_capacity=2
    )
    assert plan.total_conflicts == 0
    color = plan.assignments
    assert color[0] != color[1] != color[2] and color[0] != color[2]
    assert color[3] == color[0]  # 3 只与 1 冲突，与 0 同池新增冲突 0


def test_dsatur_position_is_final_tiebreak():
    # 无边图：全部按位置顺序，且始终放入池号最小的空池
    plan = build_pool_plan(
        primer_count=3, edges=[], pool_count=3, pool_capacity=1
    )
    assert _assignments(plan) == [0, 1, 2]
    assert plan.total_conflicts == 0


def test_dsatur_placement_prefiles_empty_pool_by_pool_number():
    # K3 三顶点三池：第三个顶点见 2 个邻池，选剩余池号；首顶点落池 0
    plan = build_pool_plan(
        primer_count=3,
        edges=[(0, 1, "high"), (1, 2, "high"), (0, 2, "high")],
        pool_count=3,
        pool_capacity=1,
    )
    assert plan.total_conflicts == 0
    assert sorted(_assignments(plan)) == [0, 1, 2]


def test_dsatur_deterministic_with_shuffled_edge_order():
    edges = [(2, 0, "high"), (0, 1, "medium"), (1, 2, "low")]
    a = build_pool_plan(primer_count=3, edges=edges, pool_count=2, pool_capacity=2)
    b = build_pool_plan(
        primer_count=3,
        edges=list(reversed(edges)),
        pool_count=2,
        pool_capacity=2,
    )
    assert a.assignments == b.assignments
    assert [(e.position_a, e.position_b, e.risk_level) for _, e in a.conflicts] == [
        (e.position_a, e.position_b, e.risk_level) for _, e in b.conflicts
    ]


def test_dsatur_conflicts_sorted_by_pool_then_positions():
    # K4 + 2 池容量 2：每池恰好一条冲突，按池号、位置排序
    plan = build_pool_plan(
        primer_count=4,
        edges=[(i, j, "high") for i in range(4) for j in range(i + 1, 4)],
        pool_count=2,
        pool_capacity=2,
    )
    assert plan.total_conflicts == 2
    keys = [(p, e.position_a, e.position_b) for p, e in plan.conflicts]
    assert keys == sorted(keys)
    assert [p for p, _ in plan.conflicts] == [0, 1]
