"""conflict-repair-plan 接口端到端测试与分支限界单测。"""

from __future__ import annotations

import itertools
import os
import random

import pytest

from app.repair import InfeasibleRepair, minimum_cost_removal


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("repair_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


# 各引物的风险画像（阈值 ct=4/tt=3，已用 analyze_pair 独立核算）：
#   F(pos0)=GGGGGAAAAA  自 low；R(pos1)=GGGGGTTTTT 自 low；F-R high(5,5)
#   H(pos2)=ACGTACGTAC  自 high(8,8) 自环；H 与 F/R 均 low
#   A(pos3)=ATACACGTCA  自 low
#   B(pos4)=AAGTAAGTGT  自 low；A-B high(4,4)、B-C high(3,3)
#   C(pos5)=TCTATTGACT  自 low；A-C high(4,4)
PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},
        {"name": "R", "sequence": "GGGGGTTTTT"},
        {"name": "H", "sequence": "ACGTACGTAC"},
        {"name": "A", "sequence": "ATACACGTCA"},
        {"name": "B", "sequence": "AAGTAAGTGT"},
        {"name": "C", "sequence": "TCTATTGACT"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _plan(client, batch_id, positions_costs, protected_names, levels=("high",)):
    return client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": list(levels),
            "primers": [
                {"position": pos, "removal_cost": cost}
                for pos, cost in positions_costs
            ],
            "protected_names": protected_names,
        },
    )


# ---------- 基本求解结果 ----------

def test_single_high_edge_cost_tie_picks_earliest_position(client, batch_id):
    # 仅 F-R 一条 high 边；代价相同 → 移除数量相同 → 位置序列取最靠前的 F(0)
    resp = _plan(
        client, batch_id, [(0, 1), (1, 1)], protected_names=[]
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["batch_id"] == batch_id
    assert data["risk_levels"] == ["high"]
    assert data["primer_count"] == 2
    assert data["conflict_count"] == 1
    assert [i["name"] for i in data["removed"]] == ["F"]
    assert data["total_removal_cost"] == 1
    assert data["removed_count"] == 1
    assert data["retained_count"] == 1
    assert [i["name"] for i in data["retained"]] == ["R"]
    assert data["retained"][0]["protected"] is False

    edge = data["covered_conflicts"][0]
    assert edge["name_a"] == "F"
    assert edge["name_b"] == "R"
    assert edge["position_a"] == 0
    assert edge["position_b"] == 1
    assert edge["is_self_dimer"] is False
    assert edge["risk_level"] == "high"
    assert edge["covered_by"] == ["F"]


def test_cheaper_endpoint_is_removed(client, batch_id):
    resp = _plan(
        client, batch_id, [(0, 10), (1, 1)], protected_names=[]
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert [i["name"] for i in data["removed"]] == ["R"]
    assert data["total_removal_cost"] == 1


def test_triangle_needs_two_removals_min_cost_pair(client, batch_id):
    # 三角形 A-B-C（A-B/A-C/B-C 均 high）：顶点覆盖至少两点，取代价最小的两个
    resp = _plan(
        client,
        batch_id,
        [(3, 5), (4, 1), (5, 1)],
        protected_names=[],
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["conflict_count"] == 3
    assert [i["name"] for i in data["removed"]] == ["B", "C"]
    assert data["total_removal_cost"] == 2
    assert data["removed_count"] == 2
    # 被覆盖边的 covered_by 均可由移除项解释
    removed_names = {i["name"] for i in data["removed"]}
    for edge in data["covered_conflicts"]:
        assert edge["covered_by"]
        assert set(edge["covered_by"]) <= removed_names


def test_triangle_cost_weighting(client, batch_id):
    # 抬高 B、C 代价后，A+B 与 A+C 同代价时按位置取 A、B
    resp = _plan(
        client,
        batch_id,
        [(3, 6), (4, 10), (5, 10)],
        protected_names=[],
    )
    data = resp.json()
    assert [i["name"] for i in data["removed"]] == ["A", "B"]
    assert data["total_removal_cost"] == 16


def test_cost_tie_count_then_position_on_triangle(client, batch_id):
    # 等代价三角形：两点覆盖，数量恒为 2，位置序列取字典序最小的 A、B
    resp = _plan(
        client,
        batch_id,
        [(3, 5), (4, 5), (5, 5)],
        protected_names=[],
    )
    data = resp.json()
    assert [i["name"] for i in data["removed"]] == ["A", "B"]


# ---------- 自环与强制移除 ----------

def test_self_loop_forces_removal(client, batch_id):
    # H 自二聚体 high 自环，即使 H 代价很高也必须移除；F-R 另一条边
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1), (2, 100)],
        protected_names=[],
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    removed_names = {i["name"] for i in data["removed"]}
    assert "H" in removed_names
    # F-R 边取位置靠前的 F
    assert removed_names == {"H", "F"}
    assert data["total_removal_cost"] == 101

    loop = next(e for e in data["covered_conflicts"] if e["is_self_dimer"])
    assert loop["position_a"] == loop["position_b"] == 2
    assert loop["name_a"] == loop["name_b"] == "H"
    assert loop["risk_level"] == "high"
    assert loop["covered_by"] == ["H"]

    # 全部冲突（自环 + F-R）都被列出且覆盖
    assert data["conflict_count"] == 2
    assert len(data["covered_conflicts"]) == 2
    for edge in data["covered_conflicts"]:
        assert edge["covered_by"]
        for name in edge["covered_by"]:
            assert name in removed_names


def test_protected_self_loop_is_409(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1), (2, 1)],
        protected_names=["H"],
    )
    assert resp.status_code == 409, resp.text
    body = resp.json()
    assert "batch_id" in body
    assert body["risk_levels"] == ["high"]
    conflicts = body["unresolvable_conflicts"]
    assert len(conflicts) == 1
    c = conflicts[0]
    assert c["position_a"] == c["position_b"] == 2
    assert c["name_a"] == c["name_b"] == "H"
    assert c["is_self_dimer"] is True
    assert c["risk_level"] == "high"
    assert c["reason"] == "protected_self_loop"


# ---------- 受保护节点 ----------

def test_protected_pair_is_409(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1)],
        protected_names=["F", "R"],
    )
    assert resp.status_code == 409
    conflicts = resp.json()["unresolvable_conflicts"]
    assert len(conflicts) == 1
    c = conflicts[0]
    assert (c["name_a"], c["name_b"]) == ("F", "R")
    assert c["is_self_dimer"] is False
    assert c["reason"] == "protected_pair"


def test_protection_forces_neighbor_removal(client, batch_id):
    # 保护 F → F-R 边只能移除 R，即使 R 代价更高
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 100)],
        protected_names=["F"],
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert [i["name"] for i in data["removed"]] == ["R"]
    assert data["total_removal_cost"] == 100
    items = {i["name"]: i for i in data["retained"]}
    assert items["F"]["protected"] is True
    assert items["F"]["removal_cost"] == 1


def test_409_lists_all_unresolvable_conflicts_sorted(client, batch_id):
    # 保护 F、R（F-R 不可消除）且保护 H（H 自环不可消除）
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1), (2, 1)],
        protected_names=["F", "R", "H"],
    )
    assert resp.status_code == 409
    conflicts = resp.json()["unresolvable_conflicts"]
    # 按 (position_a, position_b) 升序：自环 H(2,2) 与 F-R(0,1)，0,1 在前
    keys = [(c["position_a"], c["position_b"]) for c in conflicts]
    assert keys == sorted(keys)
    assert keys == [(0, 1), (2, 2)]
    reasons = {
        (c["name_a"], c["name_b"]): c["reason"] for c in conflicts
    }
    assert reasons[("F", "R")] == "protected_pair"
    assert reasons[("H", "H")] == "protected_self_loop"


# ---------- 风险等级选择 ----------

def test_medium_self_loop_only_included_when_selected(client, batch_id):
    # H 自二聚体 high：仅选 medium 时 H 不产生自环；F-R 为 high 也不产生边
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1), (2, 1)],
        protected_names=["F", "R", "H"],
        levels=("medium",),
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["conflict_count"] == 0
    assert data["removed"] == []
    assert data["total_removal_cost"] == 0
    assert {i["name"] for i in data["retained"]} == {"F", "R", "H"}
    assert data["covered_conflicts"] == []


def test_risk_levels_normalized_and_deduped(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["HIGH", " high", "medium"],
            "primers": [
                {"position": 0, "removal_cost": 1},
                {"position": 1, "removal_cost": 1},
            ],
            "protected_names": [],
        },
    )
    assert resp.status_code == 200, resp.text
    assert resp.json()["risk_levels"] == ["high", "medium"]


# ---------- 确定性决胜 ----------

def test_deterministic_and_position_series_tiebreak(client, batch_id):
    # 两次请求（引物输入顺序打乱）结果一致，且移除项按批次位置输出
    body = {
        "risk_levels": ["high"],
        "primers": [
            {"position": 5, "removal_cost": 1},
            {"position": 4, "removal_cost": 1},
            {"position": 3, "removal_cost": 1},
            {"position": 0, "removal_cost": 1},
            {"position": 1, "removal_cost": 1},
        ],
        "protected_names": [],
    }
    r1 = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan", json=body
    ).json()
    r2 = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan", json=body
    ).json()
    assert r1 == r2
    removed_positions = [i["position"] for i in r1["removed"]]
    assert removed_positions == sorted(removed_positions)
    retained_positions = [i["position"] for i in r1["retained"]]
    assert retained_positions == sorted(retained_positions)
    # 三角形 A-B-C 等代价取位置最小两点 A(3)、B(4)；F-R 取 F(0)
    assert removed_positions == [0, 3, 4]


# ---------- 不修改批次 ----------

def test_batch_not_modified(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (1, 1), (2, 9), (3, 2), (4, 2), (5, 2)],
        protected_names=["A"],
    )
    assert resp.status_code == 200
    after = client.get(f"/api/batches/{batch_id}").json()
    assert after["pair_count"] == before["pair_count"]
    assert after["primer_count"] == before["primer_count"]
    assert [
        (p["pair_index"], p["name_top"], p["name_bottom"], p["risk_level"])
        for p in after["pairs"]
    ] == [
        (p["pair_index"], p["name_top"], p["name_bottom"], p["risk_level"])
        for p in before["pairs"]
    ]


# ---------- 错误处理 ----------

def test_batch_not_found_404(client):
    resp = _plan(
        client,
        "no-such-batch",
        [(0, 1), (1, 1)],
        protected_names=["x"],
    )
    assert resp.status_code == 404


def test_position_out_of_range_422(client, batch_id):
    resp = _plan(
        client, batch_id, [(0, 1), (99, 1)], protected_names=[]
    )
    assert resp.status_code == 422
    assert "99" in resp.json()["detail"]


def test_duplicate_positions_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [(0, 1), (0, 2)],
        protected_names=[],
    )
    assert resp.status_code == 422
    assert "重复" in resp.json()["detail"]


def test_protected_name_not_in_batch_422(client, batch_id):
    resp = _plan(
        client, batch_id, [(0, 1), (1, 1)], protected_names=["ghost"]
    )
    assert resp.status_code == 422
    assert "ghost" in resp.json()["detail"]


def test_protected_name_outside_request_422(client, batch_id):
    # H 是批次引物但不在请求 primers 内
    resp = _plan(
        client, batch_id, [(0, 1), (1, 1)], protected_names=["H"]
    )
    assert resp.status_code == 422
    assert "H" in resp.json()["detail"]


def test_primer_count_bounds_422(client, batch_id):
    one = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["high"],
            "primers": [{"position": 0, "removal_cost": 1}],
            "protected_names": ["F"],
        },
    )
    assert one.status_code == 422

    twenty_nine = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["high"],
            "primers": [
                {"position": i % 6, "removal_cost": 1} for i in range(29)
            ],
            "protected_names": ["F"],
        },
    )
    assert twenty_nine.status_code == 422


def test_cost_must_be_positive_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["high"],
            "primers": [
                {"position": 0, "removal_cost": 0},
                {"position": 1, "removal_cost": 1},
            ],
            "protected_names": ["F"],
        },
    )
    assert resp.status_code == 422


def test_invalid_risk_level_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["critical"],
            "primers": [
                {"position": 0, "removal_cost": 1},
                {"position": 1, "removal_cost": 1},
            ],
            "protected_names": ["F"],
        },
    )
    assert resp.status_code == 422


def test_blank_protected_name_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/conflict-repair-plan",
        json={
            "risk_levels": ["high"],
            "primers": [
                {"position": 0, "removal_cost": 1},
                {"position": 1, "removal_cost": 1},
            ],
            "protected_names": ["   "],
        },
    )
    assert resp.status_code == 422


# ---------- 分支限界精确求解单测 ----------

def _brute(edges, costs, protected):
    """穷举全部可行覆盖，按 (代价, 数量, 位置序列) 取最优；无解返回 None。"""
    n = len(costs)
    es = {(min(a, b), max(a, b)) for a, b in edges}
    impossible = {
        e
        for e in es
        if (e[0] == e[1] and e[0] in protected)
        or (e[0] != e[1] and e[0] in protected and e[1] in protected)
    }
    if impossible:
        return None
    best = None
    for r in range(n + 1):
        for comb in itertools.combinations(range(n), r):
            s = set(comb)
            if s & protected:
                continue
            if all(a in s or b in s for a, b in es):
                key = (
                    sum(costs[x] for x in s),
                    len(s),
                    tuple(sorted(s)),
                )
                if best is None or key < best[0]:
                    best = (key, tuple(sorted(s)))
    return best[1]


def test_solver_no_conflicts_keeps_everything():
    plan = minimum_cost_removal(
        [], costs=[5, 6, 7], protected=frozenset({2})
    )
    assert plan.removed == ()
    assert plan.retained == (0, 1, 2)


def test_solver_self_loop_forced_even_if_protected_neighbor_exists():
    # 自环节点非保护必移除；其保护邻居与他人的边照常求解
    plan = minimum_cost_removal(
        [(1, 1), (0, 1), (1, 2)],
        costs=[100, 1, 100],
        protected=frozenset({0}),
    )
    # 1 自环必移除（顺带覆盖 0-1、1-2），0 保留
    assert plan.removed == (1,)
    assert plan.retained == (0, 2)


def test_solver_infeasible_reasons_are_complete_and_sorted():
    try:
        minimum_cost_removal(
            [(0, 0), (2, 2), (1, 2), (0, 1)],
            costs=[1, 1, 1],
            protected=frozenset({0, 1, 2}),
        )
    except InfeasibleRepair as exc:
        # 自环 (0,0)、(2,2)，普通边 0-1、1-2 两端均保护
        assert exc.edges == [(0, 0), (0, 1), (1, 2), (2, 2)]
    else:
        raise AssertionError("expected InfeasibleRepair")


def test_solver_cost_then_count_then_position_ordering():
    # 路径 0-1-2：移除 {1}(代价3) vs {0,2}(代价1+1=2) → {0,2} 更便宜
    plan = minimum_cost_removal(
        [(0, 1), (1, 2)], costs=[1, 3, 1], protected=frozenset()
    )
    assert plan.removed == (0, 2)

    # 两条独立边 0-1、2-3，代价相同：单边各取较小端点 → {0,2}
    plan = minimum_cost_removal(
        [(0, 1), (2, 3)], costs=[1, 1, 1, 1], protected=frozenset()
    )
    assert plan.removed == (0, 2)


def test_solver_count_tiebreak_on_star_graph():
    # 星型：中心 0 连 1、2、3。移除中心 {0} 代价 3（1 个点）；
    # 移除所有叶子 {1,2,3} 代价 1+1+1=3（3 个点）。代价相同 → 数量少者优先
    plan = minimum_cost_removal(
        [(0, 1), (0, 2), (0, 3)],
        costs=[3, 1, 1, 1],
        protected=frozenset(),
    )
    assert plan.removed == (0,)

    # 中心代价 4 时叶子方案（3）更便宜，数量决胜不再生效
    plan = minimum_cost_removal(
        [(0, 1), (0, 2), (0, 3)],
        costs=[4, 1, 1, 1],
        protected=frozenset(),
    )
    assert plan.removed == (1, 2, 3)


def test_solver_matches_brute_force_randomized():
    rng = random.Random(20260919)
    for _ in range(8000):
        n = rng.randint(2, 12)
        edges = set()
        for v in range(n):
            if rng.random() < 0.15:
                edges.add((v, v))
        for u in range(n):
            for v in range(u + 1, n):
                if rng.random() < 0.28:
                    edges.add((u, v))
        costs = [rng.randint(1, 12) for _ in range(n)]
        protected = frozenset(v for v in range(n) if rng.random() < 0.2)
        try:
            plan = minimum_cost_removal(
                list(edges), costs=costs, protected=protected
            )
        except InfeasibleRepair as exc:
            assert _brute(list(edges), costs, protected) is None
            es = {(min(a, b), max(a, b)) for a, b in edges}
            expected = sorted(
                {
                    e
                    for e in es
                    if (e[0] == e[1] and e[0] in protected)
                    or (
                        e[0] != e[1]
                        and e[0] in protected
                        and e[1] in protected
                    )
                }
            )
            assert exc.edges == expected
            continue
        expected = _brute(list(edges), costs, protected)
        assert plan.removed == expected
        for u, v in edges:
            assert u in plan.removed or v in plan.removed
        assert not (set(plan.removed) & protected)


def test_solver_scale_28_is_fast():
    import time

    edges = [(i, j) for i in range(28) for j in range(i + 1, 28)]
    start = time.time()
    plan = minimum_cost_removal(
        edges, costs=[(i % 7) + 1 for i in range(28)],
        protected=frozenset(),
    )
    assert time.time() - start < 2
    # K28 的顶点覆盖须含 27 个点；结果必须真的覆盖全部边
    removed = set(plan.removed)
    assert all(u in removed or v in removed for u, v in edges)
