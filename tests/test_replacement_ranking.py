"""替换候选排名接口端到端测试与非支配分层单测。"""

from __future__ import annotations

import os

import pytest

from app.dimer import analyze_pair
from app.ranking import dominates, non_dominated_layers


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("ranking_db")
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
# 原批次：F 自二聚体 low(0,0)；R 自二聚体 high(5,5)；F-R high(run=5,term=5)

# 候选（指标用 analyze_pair 独立核算）：
#   good GGGGGAAAA : self low(0,0)，vsR high(run=4,term=4) -> (1, 4, 4)
#   keep GGGGGAAAAA: self low(0,0)，vsR high(run=5,term=5) -> (1, 5, 5)
#   twin TTTTTCCCCC: self low(0,0)，vsR high(run=5,term=5) -> (1, 5, 5)
#   bad  CCCCCGGGGG: self high(run=10,term=10)，vsR medium(run=5,term=0)
#                                                           -> (1, 10, 10)
CANDIDATES = [
    {"label": "bad", "sequence": "CCCCCGGGGG"},
    {"label": "keep", "sequence": "GGGGGAAAAA"},
    {"label": "twin", "sequence": "TTTTTCCCCC"},
    {"label": "good", "sequence": "GGGGGAAAA"},
]
EXPECTED_METRICS = {
    "good": (1, 4, 4),
    "keep": (1, 5, 5),
    "twin": (1, 5, 5),
    "bad": (1, 10, 10),
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _rank(client, batch_id, target="F", candidates=None):
    body = {"target_name": target, "candidates": candidates or CANDIDATES}
    return client.post(
        f"/api/batches/{batch_id}/replacement-ranking", json=body
    )


# ---------- 排名结果 ----------

def test_layers_ranks_and_input_order_within_layer(client, batch_id):
    data = _rank(client, batch_id).json()

    assert data["batch_id"] == batch_id
    assert data["target_name"] == "F"
    assert data["primer_count"] == 2
    assert data["candidate_count"] == 4
    assert data["continuous_threshold"] == 4
    assert data["three_prime_threshold"] == 3

    by_label = {c["label"]: c for c in data["candidates"]}

    # good 独占第 1 层；keep/twin 同层且按输入顺序（bad 在 keep 之前提交，
    # 但 bad 落到第 3 层；同层的 keep 先于 twin）；bad 第 3 层
    assert [(c["label"], c["layer"], c["rank"]) for c in data["candidates"]] == [
        ("good", 1, 1),
        ("keep", 2, 2),
        ("twin", 2, 3),
        ("bad", 3, 4),
    ]

    for label, (high, m3, mr) in EXPECTED_METRICS.items():
        cand = by_label[label]
        assert cand["high_pair_count"] == high
        assert cand["max_three_prime_run"] == m3
        assert cand["max_continuous_run"] == mr
        assert cand["sequence"] == next(
            cc["sequence"] for cc in CANDIDATES if cc["label"] == label
        )


def test_pair_details_match_existing_algorithm_and_indices(client, batch_id):
    data = _rank(client, batch_id).json()
    by_label = {c["label"]: c for c in data["candidates"]}

    # 每个候选的配对明细仅覆盖两个受影响配对：F 自二聚体（pair_index=0）
    # 与 F-R（pair_index=2），不包含 R 的自二聚体（无关配对不重算）
    keep = by_label["keep"]
    keep_detail = {
        p["pair_index"]: p
        for p in (
            keep["high_pairs"]
            + keep["max_three_prime_pairs"]
            + keep["max_continuous_pairs"]
        )
    }
    # 自二聚体 low(0,0) 未达到任何最大值，仅异二聚体 F-R 入明细
    assert sorted(keep_detail) == [2]
    keep_hetero = keep_detail[2]
    assert keep_hetero["is_self_dimer"] is False
    assert keep_hetero["name_top"] == "F"
    assert keep_hetero["name_bottom"] == "R"
    assert keep_hetero["risk_level"] == "high"

    # 与直接调用现有 analyze_pair 的结果逐项一致
    expected = analyze_pair("F", "GGGGGAAAAA", "R", "GGGGGTTTTT", 4, 3)
    assert keep_hetero["continuous_complementary_length"] == 5
    assert keep_hetero["three_prime_run_length"] == 5
    assert keep_hetero["best_shift"] == expected.best_shift
    assert keep_hetero["alignment"]["top"] == expected.top_line
    assert keep_hetero["risk_reasons"] == list(expected.risk_reasons)

    # bad 的自二聚体 high(run=10)，异二聚体仅 medium：high_pairs 只有自二聚体
    bad = by_label["bad"]
    assert [(p["name_top"], p["name_bottom"]) for p in bad["high_pairs"]] == [
        ("F", "F")
    ]
    # 最大 3′ 连续由自二聚体达到（10）
    assert [
        (p["name_top"], p["name_bottom"])
        for p in bad["max_three_prime_pairs"]
    ] == [("F", "F")]
    # 最大连续互补自二聚体(10) 与异二聚体(5) 只有自二聚体达到 10
    assert [
        p["continuous_complementary_length"]
        for p in bad["max_continuous_pairs"]
    ] == [10]

    # good 的最大值由异二聚体达到（4），自二聚体为 0 不入选
    good = by_label["good"]
    assert [
        (p["name_top"], p["name_bottom"])
        for p in good["max_three_prime_pairs"]
    ] == [("F", "R")]
    assert good["max_three_prime_pairs"][0]["three_prime_run_length"] == 4


def test_metrics_independently_recomputed_with_analyze_pair(client, batch_id):
    data = _rank(client, batch_id).json()
    other = {"R": "GGGGGTTTTT"}
    for cand in data["candidates"]:
        seq = cand["sequence"]
        pairs = [
            analyze_pair("F", seq, "F", seq, 4, 3, is_self_dimer=True),
            analyze_pair("F", seq, "R", other["R"], 4, 3),
        ]
        assert cand["high_pair_count"] == sum(
            p.risk_level == "high" for p in pairs
        )
        assert cand["max_three_prime_run"] == max(
            p.terminal_run for p in pairs
        )
        assert cand["max_continuous_run"] == max(
            p.max_complementary_run for p in pairs
        )


def test_in_memory_only_no_database_write(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    resp = _rank(client, batch_id)
    assert resp.status_code == 200
    after = client.get(f"/api/batches/{batch_id}").json()

    # 批次、引物、配对全部保持原样
    assert after["pair_count"] == before["pair_count"]
    assert [p["name"] for p in after["primers"]] == ["F", "R"]
    assert [
        (p["pair_index"], p["name_top"], p["name_bottom"]) for p in after["pairs"]
    ] == [
        (p["pair_index"], p["name_top"], p["name_bottom"])
        for p in before["pairs"]
    ]
    # 候选序列不出现在持久化引物中
    assert "CCCCCGGGGG" not in {p["sequence"] for p in after["primers"]}

    # 批次列表数量不变
    listing = client.get("/api/batches").json()
    assert all(
        b["id"] != batch_id or b["pair_count"] == 3 for b in listing["batches"]
    )


def test_only_affected_pairs_are_recomputed(client, batch_id, monkeypatch):
    # 三引物批次：替换中间引物时应产生 1 个自二聚体 + 2 个异二聚体，共 3 次
    # analyze_pair 调用/候选；原批次的其余配对不重算
    three = client.post(
        "/api/batches",
        json={
            "primers": [
                {"name": "F", "sequence": "GGGGGAAAAA"},
                {"name": "R", "sequence": "GGGGGTTTTT"},
                {"name": "X", "sequence": "AAAAAAAAAA"},
            ],
            "continuous_threshold": 4,
            "three_prime_threshold": 3,
        },
    ).json()

    import app.services as services

    calls = []
    original = services.analyze_pair

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(services, "analyze_pair", spy)
    resp = client.post(
        f"/api/batches/{three['id']}/replacement-ranking",
        json={
            "target_name": "R",
            "candidates": [
                {"label": "a", "sequence": "ACGTACGT"},
                {"label": "b", "sequence": "TTTTGGGG"},
            ],
        },
    )
    assert resp.status_code == 200
    assert len(calls) == 2 * 3
    details_by_label = {
        cand["label"]: sorted(
            {p["pair_index"] for p in cand["high_pairs"]}
            | {p["pair_index"] for p in cand["max_three_prime_pairs"]}
            | {p["pair_index"] for p in cand["max_continuous_pairs"]}
        )
        for cand in resp.json()["candidates"]
    }
    # a=ACGTACGT：仅自二聚体 high(8,8) 达最大值；两个异二聚体 low(1,1)
    # b=TTTTGGGG：自二聚体 low(0,0)，两个异二聚体 high(4,4) 达最大值
    # 自二聚体 pair_index=1（R 的位置）；F-R=3、R-X=5；F-X=4 无关不重算
    assert details_by_label == {"a": [1], "b": [3, 5]}


def test_single_primer_batch_only_self_dimer(client):
    solo = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "solo", "sequence": "AAAATTTT"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 3,
        },
    ).json()
    resp = client.post(
        f"/api/batches/{solo['id']}/replacement-ranking",
        json={
            "target_name": "solo",
            "candidates": [
                {"label": "a", "sequence": "AAAATTTT"},
                {"label": "b", "sequence": "ACGTACGT"},
            ],
        },
    )
    assert resp.status_code == 200
    for cand in resp.json()["candidates"]:
        all_idx = {
            p["pair_index"]
            for p in (
                cand["high_pairs"]
                + cand["max_three_prime_pairs"]
                + cand["max_continuous_pairs"]
            )
        }
        assert all_idx == {0}


# ---------- 错误处理 ----------

def test_batch_not_found_404(client):
    resp = client.post(
        "/api/batches/no-such-batch/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": "a", "sequence": "ACGT"},
                {"label": "b", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 404


def test_target_primer_missing_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "ghost",
            "candidates": [
                {"label": "a", "sequence": "ACGT"},
                {"label": "b", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 422
    assert "ghost" in resp.json()["detail"]


def test_duplicate_labels_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": "dup", "sequence": "ACGT"},
                {"label": "dup", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 422
    assert "重复" in resp.json()["detail"]


def test_duplicate_labels_after_trim_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": "tag ", "sequence": "ACGT"},
                {"label": " tag", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 422


def test_candidate_count_bounds_422(client, batch_id):
    one = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [{"label": "only", "sequence": "ACGT"}],
        },
    )
    assert one.status_code == 422

    many = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": f"c{i}", "sequence": "ACGT"} for i in range(51)
            ],
        },
    )
    assert many.status_code == 422


def test_invalid_iupac_candidate_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": "a", "sequence": "ACGX"},
                {"label": "b", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 422


def test_empty_label_422(client, batch_id):
    resp = client.post(
        f"/api/batches/{batch_id}/replacement-ranking",
        json={
            "target_name": "F",
            "candidates": [
                {"label": "  ", "sequence": "ACGT"},
                {"label": "b", "sequence": "TTTT"},
            ],
        },
    )
    assert resp.status_code == 422


# ---------- 非支配分层单测 ----------

def test_dominates_smaller_is_better():
    assert dominates((0, 1, 2), (1, 1, 2))       # 一维严格更优
    assert dominates((1, 1, 1), (2, 2, 2))       # 全维严格更优
    assert not dominates((1, 2, 1), (1, 1, 1))   # 有维度更差
    assert not dominates((1, 1, 1), (1, 1, 1))   # 相同不构成支配
    assert not dominates((1, 0, 2), (0, 1, 2))   # 互有胜负


def test_non_dominated_layers_order_and_stability():
    items = ["bad", "keep", "twin", "good"]
    vec = {
        "good": (1, 4, 4),
        "keep": (1, 5, 5),
        "twin": (1, 5, 5),
        "bad": (1, 10, 10),
    }
    layers = non_dominated_layers(items, lambda x: vec[x])
    assert layers == [["good"], ["keep", "twin"], ["bad"]]

    # 同向量元素互不支配、同层；层内稳定保持输入顺序
    items2 = ["twin", "keep", "good", "bad"]
    layers2 = non_dominated_layers(items2, lambda x: vec[x])
    assert layers2 == [["good"], ["twin", "keep"], ["bad"]]


def test_non_dominated_layers_all_equal_single_layer():
    items = ["a", "b", "c"]
    layers = non_dominated_layers(items, lambda _: (0, 0, 0))
    assert layers == [["a", "b", "c"]]
