"""augmentation-plan 接口端到端测试与分支限界单测。"""

from __future__ import annotations

import itertools
import os
from itertools import combinations

import pytest

from app.augmentation import maximum_independent_set


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("aug_db")
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
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _plan(client, batch_id, candidates, max_selected=10):
    return client.post(
        f"/api/batches/{batch_id}/augmentation-plan",
        json={"max_selected": max_selected, "candidates": candidates},
    )


# 候选风险画像（阈值 ct=4/tt=3，已用 analyze_pair 独立核算）：
#   self-hairpin ACGTACGTAC : 自二聚体 high(8,8)，与 F/R 均 low -> 自二聚体阻断
#   vsF-blocked TTTTGGGGAA  : 自 low，与 F high(4,4)，与 R low -> 原批次阻断
#   其余 5 条自/对 F/R 均 low，候选间冲突（异二聚体 high）：
#     A=ATACACGTCA, B=AAGTAAGTGT, C=TCTATTGACT：A-B/A-C/B-C 两两 high（三角形）
#     D=CTTAAGGGTT 与所有候选均 low（自由点）
#     E=ATTACACTCA：B-E high(5,5)、C-E high(3,3)，与 A、D 不冲突
#   最大无冲突集大小为 3，字典序最小解为 {A, D, E}
def _cand(name, seq):
    return {"name": name, "sequence": seq}


CANDIDATES = [
    _cand("self-hairpin", "ACGTACGTAC"),
    _cand("A", "ATACACGTCA"),
    _cand("B", "AAGTAAGTGT"),
    _cand("vsF-blocked", "TTTTGGGGAA"),
    _cand("C", "TCTATTGACT"),
    _cand("D", "CTTAAGGGTT"),
    _cand("E", "ATTACACTCA"),
]


def _brute_expected(candidate_names, conflict, cap):
    """在候选间冲突图（含自环）上求字典序最小的最大独立集（容量 cap）。

    conflict: dict 名称 -> set(名称)，包含自身名即自环阻断。
    """
    names = list(candidate_names)
    free = [n for n in names if n not in conflict.get(n, ())]
    best = ()
    for r in range(min(cap, len(free)), -1, -1):
        sols = [
            comb
            for comb in combinations(free, r)
            if all(b not in conflict[a] for a, b in itertools.combinations(comb, 2))
        ]
        if sols:
            best = min(sols)
            break
    return best


# ---------- 基本选择结果 ----------

def test_response_envelope_and_thresholds(client, batch_id):
    resp = _plan(client, batch_id, CANDIDATES, max_selected=10)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["batch_id"] == batch_id
    assert data["continuous_threshold"] == 4
    assert data["three_prime_threshold"] == 3
    assert data["max_selected"] == 10
    assert data["candidate_count"] == 7
    assert data["selected_count"] == 3
    assert data["rejected_count"] == 4
    # 冲突三角形中取请求顺序最靠前的 A；D 为自由点；E 与 A、D 不冲突
    assert [c["name"] for c in data["selected"]] == ["A", "D", "E"]
    assert [c["position"] for c in data["selected"]] == [1, 5, 6]


def test_selection_matches_brute_force_maximum_conflict_free_set(client, batch_id):
    # 候选间冲突（仅 high 阻断；此处候选间阻断边均为 high）
    conflict = {
        "A": {"B", "C"},
        "B": {"A", "C", "E"},
        "C": {"A", "B", "E"},
        "D": set(),
        "E": {"B", "C"},
        "self-hairpin": {"self-hairpin"},
        "vsF-blocked": {"vsF-blocked"},
    }
    names = [c["name"] for c in CANDIDATES]
    expected = _brute_expected(names, conflict, cap=10)

    data = _plan(client, batch_id, CANDIDATES, max_selected=10).json()
    assert tuple(c["name"] for c in data["selected"]) == expected


def test_self_dimer_blocked_candidate_lists_self_pair(client, batch_id):
    data = _plan(client, batch_id, CANDIDATES, max_selected=10).json()
    rej = {c["name"]: c for c in data["rejected"]}

    hair = rej["self-hairpin"]
    assert hair["reason"] == "risk"
    assert hair["limited_by_capacity"] is False
    pairs = hair["blocking_pairs"]
    assert len(pairs) == 1
    pair = pairs[0]
    assert pair["is_self_dimer"] is True
    assert pair["counterpart_kind"] == "self"
    assert pair["counterpart_name"] == "self-hairpin"
    assert pair["name_top"] == pair["name_bottom"] == "self-hairpin"
    assert pair["risk_level"] == "high"
    assert pair["continuous_complementary_length"] == 8
    assert pair["three_prime_run_length"] == 8
    assert pair["risk_reasons"]


def test_batch_blocked_candidate_lists_batch_pair_and_metrics(client, batch_id):
    data = _plan(client, batch_id, CANDIDATES, max_selected=10).json()
    rej = {c["name"]: c for c in data["rejected"]}

    blocked = rej["vsF-blocked"]
    assert blocked["reason"] == "risk"
    assert blocked["limited_by_capacity"] is False
    pairs = blocked["blocking_pairs"]
    assert len(pairs) == 1
    pair = pairs[0]
    assert pair["is_self_dimer"] is False
    assert pair["counterpart_kind"] == "batch"
    assert pair["counterpart_name"] == "F"
    assert pair["risk_level"] == "high"
    assert pair["continuous_complementary_length"] == 4
    assert pair["three_prime_run_length"] == 4


def test_conflict_with_selected_candidate_is_listed(client, batch_id):
    data = _plan(client, batch_id, CANDIDATES, max_selected=10).json()
    rej = {c["name"]: c for c in data["rejected"]}

    # B 与入选的 A、E 冲突（B-C 边存在，但 C 未入选，故不列入阻断配对）
    b = rej["B"]
    assert b["reason"] == "risk"
    assert b["limited_by_capacity"] is False
    counterparts = {
        (p["counterpart_kind"], p["counterpart_name"]): p for p in b["blocking_pairs"]
    }
    assert set(counterparts) == {
        ("candidate", "A"),
        ("candidate", "E"),
    }
    be = counterparts[("candidate", "E")]
    assert be["name_top"] == "B"
    assert be["name_bottom"] == "E"
    assert be["risk_level"] == "high"
    assert be["continuous_complementary_length"] == 5
    assert be["three_prime_run_length"] == 5

    # 明细按全局配对序号（pair_index）升序
    indices = [p["pair_index"] for p in b["blocking_pairs"]]
    assert indices == sorted(indices)


# ---------- 容量上限 ----------

def test_capacity_limit_marks_capacity_rejections(client, batch_id):
    safe_candidates = [
        _cand("A", "ATACACGTCA"),
        _cand("D", "CTTAAGGGTT"),
        _cand("E", "ATTACACTCA"),
    ]
    data = _plan(client, batch_id, safe_candidates, max_selected=2).json()

    assert data["selected_count"] == 2
    assert data["rejected_count"] == 1
    # 同规模字典序取请求顺序最靠前的两个
    assert [c["name"] for c in data["selected"]] == ["A", "D"]

    e = next(c for c in data["rejected"] if c["name"] == "E")
    assert e["reason"] == "capacity"
    assert e["limited_by_capacity"] is True
    assert e["blocking_pairs"] == []


def test_capacity_below_maximum_favors_lexicographically_earliest(client, batch_id):
    # 冲突三角形 A-B-C：容量 1 时取请求顺序最靠前且可独立存在的 A
    triangle = [
        _cand("A", "ATACACGTCA"),
        _cand("B", "AAGTAAGTGT"),
        _cand("C", "TCTATTGACT"),
    ]
    data = _plan(client, batch_id, triangle, max_selected=1).json()
    assert [c["name"] for c in data["selected"]] == ["A"]
    for name in ("B", "C"):
        rej = next(c for c in data["rejected"] if c["name"] == name)
        assert rej["reason"] == "risk"
        # 阻断配对必须指向入选候选 A
        assert all(
            p["counterpart_kind"] == "candidate" and p["counterpart_name"] == "A"
            for p in rej["blocking_pairs"]
        )


def test_reordered_request_flips_lexicographic_choice(client, batch_id):
    triangle_first = [
        _cand("B", "AAGTAAGTGT"),
        _cand("A", "ATACACGTCA"),
        _cand("C", "TCTATTGACT"),
    ]
    data = _plan(client, batch_id, triangle_first, max_selected=1).json()
    # 请求顺序变为 B 最先时，最优解（容量 1）唯一地取 B
    assert [c["name"] for c in data["selected"]] == ["B"]


def test_capacity_larger_than_feasible_does_not_limit(client, batch_id):
    data = _plan(client, batch_id, CANDIDATES, max_selected=18).json()
    assert data["selected_count"] == 3
    assert all(
        c["reason"] == "risk" for c in data["rejected"]
    )


def test_medium_only_conflict_is_also_blocking(client, batch_id):
    # M1/M2 自二聚体与对 F/R 均 low，仅相互之间为 medium（连续 4、3′端 1）
    medium_pair = [
        _cand("M1", "AATGTAGGCG"),
        _cand("M2", "ACCCCTATAA"),
    ]
    data = _plan(client, batch_id, medium_pair, max_selected=1).json()
    assert [c["name"] for c in data["selected"]] == ["M1"]

    m2 = next(c for c in data["rejected"] if c["name"] == "M2")
    assert m2["reason"] == "risk"
    assert m2["limited_by_capacity"] is False
    assert len(m2["blocking_pairs"]) == 1
    pair = m2["blocking_pairs"][0]
    assert pair["risk_level"] == "medium"
    assert pair["counterpart_kind"] == "candidate"
    assert pair["counterpart_name"] == "M1"
    assert pair["continuous_complementary_length"] == 4
    assert pair["three_prime_run_length"] == 1

    # 容量放宽到 2 时仍无法同时入选（medium 同样构成冲突边）
    data2 = _plan(client, batch_id, medium_pair, max_selected=2).json()
    assert [c["name"] for c in data2["selected"]] == ["M1"]


def test_multiple_blocking_pairs_self_and_batch_medium(client, batch_id):
    # mild=CCCCGGGGTT：自二聚体 medium(8,0)，与 F low，与 R medium(4,0)
    # D=CTTAAGGGTT：自/对批次均 low，与 mild 不冲突，可正常入选
    data = _plan(
        client,
        batch_id,
        [
            _cand("mild", "CCCCGGGGTT"),
            _cand("D", "CTTAAGGGTT"),
        ],
        max_selected=2,
    ).json()

    assert [c["name"] for c in data["selected"]] == ["D"]
    mild = next(c for c in data["rejected"] if c["name"] == "mild")
    assert mild["reason"] == "risk"
    pairs = mild["blocking_pairs"]
    # 自二聚体在前、原批次配对按批次位置随后
    assert [p["counterpart_kind"] for p in pairs] == ["self", "batch"]
    assert pairs[0]["is_self_dimer"] is True
    assert pairs[0]["risk_level"] == "medium"
    assert pairs[0]["continuous_complementary_length"] == 8
    assert pairs[1]["counterpart_name"] == "R"
    assert pairs[1]["risk_level"] == "medium"


# ---------- 不落库 ----------

def test_in_memory_only_no_database_write(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    resp = _plan(client, batch_id, CANDIDATES, max_selected=10)
    assert resp.status_code == 200
    after = client.get(f"/api/batches/{batch_id}").json()

    assert after["pair_count"] == before["pair_count"]
    assert [p["name"] for p in after["primers"]] == ["F", "R"]
    assert [
        (p["pair_index"], p["name_top"], p["name_bottom"]) for p in after["pairs"]
    ] == [
        (p["pair_index"], p["name_top"], p["name_bottom"]) for p in before["pairs"]
    ]

    listing = client.get("/api/batches").json()
    assert all(b["id"] != batch_id or b["pair_count"] == 3 for b in listing["batches"])


# ---------- 错误处理 ----------

def test_batch_not_found_404(client):
    resp = _plan(
        client,
        "no-such-batch",
        [_cand("a", "ACGT"), _cand("b", "TTTT")],
        max_selected=2,
    )
    assert resp.status_code == 404


def test_duplicate_candidate_names_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_cand("dup", "ACGT"), _cand("dup", "TTTT")],
        max_selected=2,
    )
    assert resp.status_code == 422
    assert "重复" in resp.json()["detail"]


def test_duplicate_names_after_trim_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_cand(" tag", "ACGT"), _cand("tag ", "TTTT")],
        max_selected=2,
    )
    assert resp.status_code == 422


def test_candidate_name_clashing_with_batch_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_cand("F", "ACGTACGT"), _cand("ok", "TTTTAAAA")],
        max_selected=2,
    )
    assert resp.status_code == 422
    assert "F" in resp.json()["detail"]


def test_candidate_count_bounds_422(client, batch_id):
    one = _plan(
        client, batch_id, [_cand("only", "ACGT")], max_selected=1
    )
    assert one.status_code == 422

    nineteen = _plan(
        client,
        batch_id,
        [_cand(f"c{i}", "ACGT") for i in range(19)],
        max_selected=18,
    )
    assert nineteen.status_code == 422


def test_max_selected_bounds_422(client, batch_id):
    cands = [_cand("a", "ACGT"), _cand("b", "TTTT")]
    assert _plan(client, batch_id, cands, max_selected=0).status_code == 422
    assert _plan(client, batch_id, cands, max_selected=19).status_code == 422


def test_invalid_iupac_candidate_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_cand("a", "ACGX"), _cand("b", "TTTT")],
        max_selected=2,
    )
    assert resp.status_code == 422


def test_empty_name_422(client, batch_id):
    resp = _plan(
        client,
        batch_id,
        [_cand("   ", "ACGT"), _cand("b", "TTTT")],
        max_selected=2,
    )
    assert resp.status_code == 422


# ---------- 分支限界精确求解单测 ----------

def _brute_adj(adj, cap):
    n = len(adj)
    verts = [v for v in range(n) if v not in adj[v]]
    for r in range(min(cap, len(verts)), -1, -1):
        sols = [
            comb
            for comb in itertools.combinations(verts, r)
            if all(b not in adj[a] for a, b in itertools.combinations(comb, 2))
        ]
        if sols:
            return min(sols)
    return ()


def test_mis_empty_graph_respects_capacity():
    adj = [set() for _ in range(5)]
    assert maximum_independent_set(adj, 3).selected == (0, 1, 2)
    assert maximum_independent_set(adj, 5).selected == (0, 1, 2, 3, 4)


def test_mis_self_loops_are_excluded():
    adj = [set() for _ in range(4)]
    adj[1].add(1)
    adj[3].add(3)
    res = maximum_independent_set(adj, 10)
    assert res.selected == (0, 2)
    assert res.self_blocked == frozenset({1, 3})


def test_mis_triangle_picks_lexicographically_smallest():
    adj = [set() for _ in range(3)]
    for a, b in itertools.combinations(range(3), 2):
        adj[a].add(b)
        adj[b].add(a)
    assert maximum_independent_set(adj, 10).selected == (0,)
    assert maximum_independent_set(adj, 1).selected == (0,)


def test_mis_matches_brute_force_on_fixed_graphs():
    # 彼得森图（α=4）+ 容量裁剪
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
        (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
    ]
    adj = [set() for _ in range(10)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    for cap in (1, 2, 3, 4, 5):
        assert maximum_independent_set(adj, cap).selected == _brute_adj(adj, cap)


def test_mis_matches_brute_force_randomized():
    import random

    rng = random.Random(20260919)
    for _ in range(3000):
        n = rng.randint(0, 18)
        adj = [set() for _ in range(n)]
        for i in range(n):
            if rng.random() < 0.2:
                adj[i].add(i)
            for j in range(i + 1, n):
                if rng.random() < rng.random():
                    adj[i].add(j)
                    adj[j].add(i)
        cap = rng.randint(1, max(1, n))
        res = maximum_independent_set(adj, cap)
        assert res.selected == _brute_adj(adj, cap)
        # 自环节点永不入选
        assert not (set(res.selected) & res.self_blocked)
