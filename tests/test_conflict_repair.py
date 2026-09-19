"""conflict-repair-plan 接口端到端测试与最小代价顶点覆盖求解器单测。"""

from __future__ import annotations

import itertools
import os

import pytest

from app.repair import minimum_cost_cover


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("repair_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


# 批次风险画像（阈值 ct=4/tt=3，已用 analyze_pair 独立核算）：
#   自二聚体：HP high(8,8)；F/R/D/M1/M2 均 low
#   异二聚体：F-R high(5,5)、HP-M1 high(3,3)、D-M2 high(3,3)、M1-M2 medium(4,1)，
#             其余两两组合均 low
# 冲突图（high+medium）：自环 HP；边 F-R、HP-M1、D-M2、M1-M2
# 冲突图（仅 high）    ：自环 HP；边 F-R、HP-M1、D-M2
PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},    # pos 0
        {"name": "R", "sequence": "GGGGGTTTTT"},    # pos 1
        {"name": "HP", "sequence": "ACGTACGTAC"},   # pos 2
        {"name": "D", "sequence": "CTTAAGGGTT"},    # pos 3
        {"name": "M1", "sequence": "AATGTAGGCG"},   # pos 4
        {"name": "M2", "sequence": "ACCCCTATAA"},   # pos 5
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _primer(name, cost=1):
    return {"name": name, "removal_cost": cost}


def _plan(client, batch_id, primers, risk_levels=None, protected=None):
    body = {
        "risk_levels": risk_levels if risk_levels is not None else ["high", "medium"],
        "primers": primers,
    }
    if protected is not None:
        body["protected_names"] = protected
    return client.post(f"/api/batches/{batch_id}/conflict-repair-plan", json=body)


# ---------- 基本求解结果 ----------

def test_response_envelope_and_min_cost_cover(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [
            _primer("F", 5), _primer("R", 1), _primer("HP", 2),
            _primer("D", 1), _primer("M1", 1), _primer("M2", 4),
        ],
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["batch_id"] == batch_id
    assert data["risk_levels"] == ["high", "medium"]
    assert data["primer_count"] == 6
    assert data["removed_count"] == 4
    assert data["kept_count"] == 2
    assert data["total_cost"] == 5

    # HP 自环强制；F-R 取便宜的 R；HP-M1 已被 HP 覆盖；
    # M1-M2 取 M1(1) 而非 M2(4)；D-M2 取 D(1) 而非 M2(4)
    assert [(p["position"], p["name"]) for p in data["removed"]] == [
        (1, "R"), (2, "HP"), (3, "D"), (4, "M1"),
    ]
    reasons = {p["name"]: p["reason"] for p in data["removed"]}
    assert reasons["HP"] == "self_loop"
    assert reasons["R"] == reasons["D"] == reasons["M1"] == "branch_and_bound"

    assert [(p["position"], p["name"], p["protected"]) for p in data["kept"]] == [
        (0, "F", False), (5, "M2", False),
    ]
    # 移除项代价之和恒等于 total_cost
    assert sum(p["removal_cost"] for p in data["removed"]) == data["total_cost"]


def test_covered_conflicts_are_annotated(client, batch_id):
    data = _plan(
        client,
        batch_id,
        [
            _primer("F", 5), _primer("R", 1), _primer("HP", 2),
            _primer("D", 1), _primer("M1", 1), _primer("M2", 4),
        ],
    ).json()

    conflicts = data["covered_conflicts"]
    # 全部 5 条冲突，按两端批次位置升序
    assert [(c["kind"], c["position_a"], c["position_b"]) for c in conflicts] == [
        ("edge", 0, 1),
        ("self_loop", 2, 2),
        ("edge", 2, 4),
        ("edge", 3, 5),
        ("edge", 4, 5),
    ]
    by_pair = {(c["name_a"], c["name_b"]): c for c in conflicts}
    assert by_pair[("F", "R")]["risk_level"] == "high"
    assert by_pair[("F", "R")]["covered_by_names"] == ["R"]
    assert by_pair[("HP", "HP")]["covered_by_names"] == ["HP"]
    assert by_pair[("HP", "M1")]["covered_by_names"] == ["HP", "M1"]
    assert by_pair[("D", "M2")]["covered_by_names"] == ["D"]
    medium = by_pair[("M1", "M2")]
    assert medium["risk_level"] == "medium"
    assert medium["covered_by_names"] == ["M1"]
    # 每条冲突至少被一个移除端点覆盖
    assert all(c["covered_by_positions"] for c in conflicts)


def test_cost_dominates_position(client, batch_id):
    # F 便宜但位置靠前、R 贵：只比成本，移除 F
    data = _plan(client, batch_id, [_primer("F", 1), _primer("R", 100)]).json()
    assert [p["name"] for p in data["removed"]] == ["F"]
    assert data["total_cost"] == 1


def test_count_tie_break_prefers_fewer_removals(client, batch_id):
    # 路径 HP-M1-M2-D（HP 自环强制后剩星形 M1-M2-D）：
    # {M2} 成本 6、1 个顶点 vs {M1,D} 成本 6、2 个顶点 -> 同成本取数量更少者
    data = _plan(
        client,
        batch_id,
        [_primer("HP", 1), _primer("D", 3), _primer("M1", 3), _primer("M2", 6)],
    ).json()
    assert [p["name"] for p in data["removed"]] == ["HP", "M2"]
    assert data["total_cost"] == 7


def test_position_sequence_tie_break(client, batch_id):
    # 两条不相交边 F-R、D-M2，全部等代价：逐边取批次位置更小端
    data = _plan(
        client,
        batch_id,
        [_primer("F"), _primer("R"), _primer("D"), _primer("M2")],
        risk_levels=["high"],
    ).json()
    assert [(p["position"], p["name"]) for p in data["removed"]] == [
        (0, "F"), (3, "D"),
    ]


def test_request_order_does_not_change_plan(client, batch_id):
    forward = _plan(
        client, batch_id, [_primer("F"), _primer("R"), _primer("D"), _primer("M2")],
        risk_levels=["high"],
    ).json()
    backward = _plan(
        client, batch_id, [_primer("M2"), _primer("D"), _primer("R"), _primer("F")],
        risk_levels=["high"],
    ).json()
    # 决胜按批次位置序列，与请求顺序无关
    assert forward["removed"] == backward["removed"]
    assert forward["covered_conflicts"] == backward["covered_conflicts"]


# ---------- 自环与受保护约束 ----------

def test_self_loop_forces_removal(client, batch_id):
    data = _plan(client, batch_id, [_primer("HP", 9), _primer("D")]).json()
    # HP 自二聚体 high：无论代价多高都必须移除
    assert [p["name"] for p in data["removed"]] == ["HP"]
    assert data["removed"][0]["reason"] == "self_loop"
    self_conflict = next(
        c for c in data["covered_conflicts"] if c["kind"] == "self_loop"
    )
    assert self_conflict["position_a"] == self_conflict["position_b"] == 2
    assert self_conflict["covered_by_names"] == ["HP"]


def test_protected_primer_is_kept_and_forces_neighbor(client, batch_id):
    data = _plan(
        client,
        batch_id,
        [_primer("F", 3), _primer("R", 1)],
        protected=["R"],
    ).json()
    # R 受保护：F-R 边只能由 F 覆盖（尽管 R 更便宜）
    assert [(p["name"], p["reason"]) for p in data["removed"]] == [
        ("F", "protected_neighbor"),
    ]
    assert data["kept"] == [
        {"position": 1, "name": "R", "removal_cost": 1, "protected": True}
    ]


def test_409_when_both_edge_endpoints_protected(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_primer("F"), _primer("R"), _primer("HP")],
        protected=["F", "R"],
    )
    assert resp.status_code == 409
    detail = resp.json()["detail"]
    assert detail["unresolvable_conflicts"] == [
        {
            "kind": "edge",
            "position_a": 0,
            "position_b": 1,
            "name_a": "F",
            "name_b": "R",
            "risk_level": "high",
        }
    ]


def test_409_when_self_loop_protected(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_primer("HP"), _primer("D")],
        protected=["HP"],
    )
    assert resp.status_code == 409
    conflicts = resp.json()["detail"]["unresolvable_conflicts"]
    assert conflicts == [
        {
            "kind": "self_loop",
            "position_a": 2,
            "position_b": 2,
            "name_a": "HP",
            "name_b": "HP",
            "risk_level": "high",
        }
    ]


def test_409_lists_all_unresolvable_conflicts(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_primer("F"), _primer("R"), _primer("HP"), _primer("M1")],
        protected=["F", "R", "HP", "M1"],
    )
    assert resp.status_code == 409
    conflicts = resp.json()["detail"]["unresolvable_conflicts"]
    # F-R 边两端受保护 + HP 自环受保护 + HP-M1 边两端受保护，按位置升序
    assert [(c["kind"], c["name_a"], c["name_b"]) for c in conflicts] == [
        ("edge", "F", "R"),
        ("self_loop", "HP", "HP"),
        ("edge", "HP", "M1"),
    ]


def test_protected_primer_without_conflicts_is_simply_kept(client, batch_id):
    data = _plan(
        client,
        batch_id,
        [_primer("D"), _primer("M1")],
        protected=["D", "M1"],
    ).json()
    # D-M1 为 low，不产生冲突：全部保留、无移除、无可覆盖冲突
    assert data["removed"] == []
    assert data["total_cost"] == 0
    assert data["covered_conflicts"] == []
    assert all(p["protected"] for p in data["kept"])


# ---------- 风险等级选择与冲突图构造 ----------

def test_risk_level_selection_changes_graph(client, batch_id):
    # M1-M2 仅为 medium：只选 high 时无冲突，选 medium 时出现冲突边
    high_only = _plan(
        client, batch_id, [_primer("M1", 1), _primer("M2", 2)],
        risk_levels=["high"],
    ).json()
    assert high_only["removed"] == []
    assert high_only["covered_conflicts"] == []

    with_medium = _plan(
        client, batch_id, [_primer("M1", 1), _primer("M2", 2)],
        risk_levels=["medium"],
    ).json()
    assert [p["name"] for p in with_medium["removed"]] == ["M1"]
    assert with_medium["total_cost"] == 1
    assert with_medium["covered_conflicts"][0]["risk_level"] == "medium"


def test_low_risk_level_can_form_edges(client, batch_id):
    # F-HP 为 low：选 low 时构成冲突边
    data = _plan(
        client, batch_id, [_primer("F", 1), _primer("HP", 2)],
        risk_levels=["low"],
    ).json()
    assert [p["name"] for p in data["removed"]] == ["F"]
    assert data["covered_conflicts"][0]["risk_level"] == "low"


def test_primer_subset_ignores_outside_pairs(client, batch_id):
    # 只指定 F、R：HP 的自环及其余配对不进入冲突图
    data = _plan(client, batch_id, [_primer("F", 2), _primer("R", 1)]).json()
    assert data["primer_count"] == 2
    assert [p["name"] for p in data["removed"]] == ["R"]
    assert [(c["name_a"], c["name_b"]) for c in data["covered_conflicts"]] == [
        ("F", "R"),
    ]


def test_risk_levels_are_normalized(client, batch_id):
    data = _plan(
        client, batch_id, [_primer("M1", 1), _primer("M2", 2)],
        risk_levels=[" Medium ", "medium", "HIGH"],
    ).json()
    # 大小写与空白规范化、重复值去重
    assert data["risk_levels"] == ["medium", "high"]
    assert [p["name"] for p in data["removed"]] == ["M1"]


# ---------- 不落库 ----------

def test_in_memory_only_no_database_write(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    resp = _plan(
        client,
        batch_id,
        [
            _primer("F", 5), _primer("R", 1), _primer("HP", 2),
            _primer("D", 1), _primer("M1", 1), _primer("M2", 4),
        ],
    )
    assert resp.status_code == 200
    after = client.get(f"/api/batches/{batch_id}").json()

    assert after["pair_count"] == before["pair_count"]
    assert [p["name"] for p in after["primers"]] == [
        p["name"] for p in before["primers"]
    ]
    assert [
        (p["pair_index"], p["name_top"], p["name_bottom"]) for p in after["pairs"]
    ] == [
        (p["pair_index"], p["name_top"], p["name_bottom"]) for p in before["pairs"]
    ]

    listing = client.get("/api/batches").json()
    assert all(
        b["id"] != batch_id or b["pair_count"] == before["pair_count"]
        for b in listing["batches"]
    )


# ---------- 错误处理 ----------

def test_batch_not_found_404(client):
    resp = _plan(client, "no-such-batch", [_primer("a"), _primer("b")])
    assert resp.status_code == 404


def test_unknown_primer_name_422(client, batch_id):
    resp = _plan(client, batch_id, [_primer("F"), _primer("ghost")])
    assert resp.status_code == 422
    assert "ghost" in resp.json()["detail"]


def test_duplicate_primer_names_422(client, batch_id):
    resp = _plan(client, batch_id, [_primer("F"), _primer(" F ")])
    assert resp.status_code == 422
    assert "重复" in resp.json()["detail"]


def test_primer_count_bounds_422(client, batch_id):
    one = _plan(client, batch_id, [_primer("F")])
    assert one.status_code == 422

    twenty_nine = _plan(
        client, batch_id, [_primer(f"p{i}") for i in range(29)]
    )
    assert twenty_nine.status_code == 422


def test_removal_cost_must_be_positive_422(client, batch_id):
    zero = _plan(client, batch_id, [_primer("F", 0), _primer("R")])
    assert zero.status_code == 422
    negative = _plan(client, batch_id, [_primer("F", -3), _primer("R")])
    assert negative.status_code == 422


def test_invalid_risk_level_422(client, batch_id):
    resp = _plan(
        client, batch_id, [_primer("F"), _primer("R")], risk_levels=["critical"]
    )
    assert resp.status_code == 422

    empty = _plan(client, batch_id, [_primer("F"), _primer("R")], risk_levels=[])
    assert empty.status_code == 422


def test_protected_name_not_in_primers_422(client, batch_id):
    resp = _plan(
        client, batch_id, [_primer("F"), _primer("R")], protected=["HP"]
    )
    assert resp.status_code == 422
    assert "HP" in resp.text


def test_duplicate_protected_names_422(client, batch_id):
    resp = _plan(
        client, batch_id, [_primer("F"), _primer("R")],
        protected=["F", " F "],
    )
    assert resp.status_code == 422
    assert "重复" in resp.text


def test_empty_names_422(client, batch_id):
    blank_primer = _plan(client, batch_id, [_primer("   "), _primer("R")])
    assert blank_primer.status_code == 422

    blank_protected = _plan(
        client, batch_id, [_primer("F"), _primer("R")], protected=["  "]
    )
    assert blank_protected.status_code == 422


# ---------- 求解器单测 ----------

def _brute_cover(costs, protected, self_loops, edges):
    """暴力枚举所有移除集，返回最优 (成本, 数量, 顶点元组) 对应的顶点元组。"""
    n = len(costs)
    if self_loops & protected:
        return None
    if any(u in protected and v in protected for u, v in edges):
        return None
    best_key = None
    best = None
    free = [v for v in range(n) if v not in protected]
    for r in range(len(free) + 1):
        for comb in itertools.combinations(free, r):
            s = set(comb)
            if not self_loops <= s:
                continue
            if not all(u in s or v in s for u, v in edges):
                continue
            key = (sum(costs[v] for v in s), len(s), tuple(sorted(s)))
            if best_key is None or key < best_key:
                best_key = key
                best = tuple(sorted(s))
    return best


def test_mcc_empty_graph_removes_nothing():
    vc = minimum_cost_cover(
        costs=[3, 1, 2], protected=frozenset(),
        self_loops=frozenset(), edges=[],
    )
    assert vc.removed == ()
    assert vc.forced == frozenset()


def test_mcc_self_loops_are_forced():
    vc = minimum_cost_cover(
        costs=[1, 1, 1], protected=frozenset(),
        self_loops=frozenset({0, 2}), edges=[(0, 1)],
    )
    assert vc.removed == (0, 2)
    assert vc.forced == frozenset({0, 2})


def test_mcc_protected_neighbor_is_forced():
    vc = minimum_cost_cover(
        costs=[1, 5], protected=frozenset({1}),
        self_loops=frozenset(), edges=[(0, 1)],
    )
    assert vc.removed == (0,)
    assert vc.forced == frozenset({0})


def test_mcc_infeasible_cases_return_none():
    assert minimum_cost_cover(
        costs=[1], protected=frozenset({0}),
        self_loops=frozenset({0}), edges=[],
    ) is None
    assert minimum_cost_cover(
        costs=[1, 1], protected=frozenset({0, 1}),
        self_loops=frozenset(), edges=[(0, 1)],
    ) is None


def test_mcc_triangle_picks_cheapest_pair_then_lexicographic():
    edges = [(0, 1), (0, 2), (1, 2)]
    vc = minimum_cost_cover(
        costs=[1, 1, 1], protected=frozenset(),
        self_loops=frozenset(), edges=edges,
    )
    assert vc.removed == (0, 1)
    vc = minimum_cost_cover(
        costs=[9, 1, 1], protected=frozenset(),
        self_loops=frozenset(), edges=edges,
    )
    assert vc.removed == (1, 2)


def test_mcc_count_tie_break_on_star():
    # 星形 K1,2：中心 6 vs 两叶 3+3，同成本取数量更少
    vc = minimum_cost_cover(
        costs=[6, 3, 3], protected=frozenset(),
        self_loops=frozenset(), edges=[(0, 1), (0, 2)],
    )
    assert vc.removed == (0,)
    # 中心 7 > 6：总成本优先，取两片叶子
    vc = minimum_cost_cover(
        costs=[7, 3, 3], protected=frozenset(),
        self_loops=frozenset(), edges=[(0, 1), (0, 2)],
    )
    assert vc.removed == (1, 2)


def test_mcc_matches_brute_force_randomized():
    import random

    rng = random.Random(20260919)
    for _ in range(3000):
        n = rng.randint(2, 11)
        costs = [rng.randint(1, 9) for _ in range(n)]
        protected = frozenset(v for v in range(n) if rng.random() < 0.2)
        self_loops = frozenset(v for v in range(n) if rng.random() < 0.15)
        edges = [
            (i, j)
            for i in range(n)
            for j in range(i + 1, n)
            if rng.random() < 0.3
        ]
        vc = minimum_cost_cover(
            costs=costs, protected=protected,
            self_loops=self_loops, edges=edges,
        )
        expected = _brute_cover(costs, protected, self_loops, edges)
        got = None if vc is None else vc.removed
        assert got == expected
        if vc is not None:
            # 自环节点必被移除，受保护节点永不被移除
            assert self_loops <= set(vc.removed)
            assert not (set(vc.removed) & protected)
            assert vc.forced <= set(vc.removed)
