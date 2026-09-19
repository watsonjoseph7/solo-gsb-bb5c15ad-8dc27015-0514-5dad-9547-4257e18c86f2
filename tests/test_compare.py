"""批次对比接口端到端测试与配对键/分类逻辑单测。"""

from __future__ import annotations

import os

import pytest

from app.compare import compare_pair_risks, pair_key


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("compare_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


THRESHOLDS = {"continuous_threshold": 4, "three_prime_threshold": 3}

# base → target 各共有配对的风险变化（已经 dimer 算法核实）：
#   A1|A2: medium→high（escalated，A2 序列变化）
#   A1|D2: high→low、D1|D2: high→low（deescalated，D2 序列变化）
#   其余共有配对等级不变；Gone 仅 base、New 仅 target
BASE_PAYLOAD = {
    "primers": [
        {"name": "A1", "sequence": "GGGGAAAAGGGG"},
        {"name": "A2", "sequence": "AAAATTTTAAAA"},
        {"name": "D1", "sequence": "GGGGGAAAAA"},
        {"name": "D2", "sequence": "GGGGGTTTTT"},
        {"name": "Gone", "sequence": "TTTTCCCC"},
    ],
    **THRESHOLDS,
}
TARGET_PAYLOAD = {
    "primers": [
        {"name": "A1", "sequence": "GGGGAAAAGGGG"},
        {"name": "A2", "sequence": "AAAATTTT"},       # 序列变化
        {"name": "D1", "sequence": "GGGGGAAAAA"},
        {"name": "D2", "sequence": "GGGGGAAAAA"},     # 序列变化
        {"name": "New", "sequence": "CCCCCCCC"},
    ],
    **THRESHOLDS,
}


@pytest.fixture(scope="module")
def batch_ids(client):
    base = client.post("/api/batches", json=BASE_PAYLOAD).json()["id"]
    target = client.post("/api/batches", json=TARGET_PAYLOAD).json()["id"]
    return base, target


def _compare(client, base_id, target_id):
    return client.get(f"/api/batches/{base_id}/compare/{target_id}")


# ---------- 四类明细与汇总 ----------

def test_compare_summary_counts(client, batch_ids):
    base_id, target_id = batch_ids
    resp = _compare(client, base_id, target_id)
    assert resp.status_code == 200, resp.text
    summary = resp.json()["summary"]
    assert summary == {
        "base_batch_id": base_id,
        "target_batch_id": target_id,
        "base_pair_count": 15,      # 5 自二聚体 + 10 异二聚体
        "target_pair_count": 15,
        "common_pair_count": 10,    # A1/A2/D1/D2 的 4 自 + 6 异
        "unchanged_pair_count": 7,
        "escalated_count": 1,
        "deescalated_count": 2,
        "introduced_count": 5,
        "resolved_count": 5,
    }


def test_compare_escalated_and_deescalated_details(client, batch_ids):
    data = _compare(client, *batch_ids).json()

    assert [p["pair_key"] for p in data["escalated"]] == ["A1|A2"]
    esc = data["escalated"][0]
    assert esc["base_risk_level"] == "medium"
    assert esc["target_risk_level"] == "high"
    assert esc["is_self_dimer"] is False
    assert esc["sequence_changed"] is True   # A2 序列变化，仍按同一配对比较

    # 按配对名确定性排序
    assert [p["pair_key"] for p in data["deescalated"]] == ["A1|D2", "D1|D2"]
    for de in data["deescalated"]:
        assert de["base_risk_level"] == "high"
        assert de["target_risk_level"] == "low"
        assert de["sequence_changed"] is True  # D2 序列变化


def test_compare_introduced_and_resolved_details(client, batch_ids):
    data = _compare(client, *batch_ids).json()

    # 自二聚体保留单名，异二聚体忽略提交顺序按名称排序
    assert [p["pair_key"] for p in data["introduced"]] == [
        "A1|New", "A2|New", "D1|New", "D2|New", "New",
    ]
    assert [p["pair_key"] for p in data["resolved"]] == [
        "A1|Gone", "A2|Gone", "D1|Gone", "D2|Gone", "Gone",
    ]

    for p in data["introduced"]:
        assert p["base_risk_level"] is None
        assert p["target_risk_level"] in {"high", "medium", "low"}
    for p in data["resolved"]:
        assert p["target_risk_level"] is None
        assert p["base_risk_level"] in {"high", "medium", "low"}

    # 明细中的自二聚体/异二聚体标记
    new_self = data["introduced"][-1]
    assert new_self["is_self_dimer"] is True
    assert new_self["name_a"] == new_self["name_b"] == "New"
    assert all(not p["is_self_dimer"] for p in data["introduced"][:-1])

    # sequence_changed 只取决于配对中同名引物的序列是否变化：
    # A2/D2 序列已改，A1/D1/New/Gone 未变
    changed = {p["pair_key"] for p in data["introduced"] if p["sequence_changed"]}
    assert changed == {"A2|New", "D2|New"}
    changed = {p["pair_key"] for p in data["resolved"] if p["sequence_changed"]}
    assert changed == {"A2|Gone", "D2|Gone"}


def test_compare_category_lists_sorted_by_pair_name(client, batch_ids):
    data = _compare(client, *batch_ids).json()
    for category in ("escalated", "deescalated", "introduced", "resolved"):
        keys = [(p["name_a"], p["name_b"]) for p in data[category]]
        assert keys == sorted(keys)


def test_compare_is_deterministic(client, batch_ids):
    a = _compare(client, *batch_ids).json()
    b = _compare(client, *batch_ids).json()
    assert a == b


# ---------- 配对键稳定性 ----------

def test_compare_ignores_submission_order(client):
    # 相同引物、不同提交顺序：全部配对键相同，无任何变化
    fwd = client.post(
        "/api/batches",
        json={
            "primers": [
                {"name": "P", "sequence": "GGGGGAAAAA"},
                {"name": "Q", "sequence": "GGGGGTTTTT"},
            ],
            **THRESHOLDS,
        },
    ).json()["id"]
    rev = client.post(
        "/api/batches",
        json={
            "primers": [
                {"name": "Q", "sequence": "GGGGGTTTTT"},
                {"name": "P", "sequence": "GGGGGAAAAA"},
            ],
            **THRESHOLDS,
        },
    ).json()["id"]

    data = _compare(client, fwd, rev).json()
    summary = data["summary"]
    assert summary["common_pair_count"] == 3     # P、Q 自二聚体 + P|Q
    assert summary["unchanged_pair_count"] == 3
    for category in ("escalated", "deescalated", "introduced", "resolved"):
        assert data[category] == []


def test_compare_same_batch_is_all_unchanged(client, batch_ids):
    base_id, _ = batch_ids
    data = _compare(client, base_id, base_id).json()
    assert data["summary"]["common_pair_count"] == 15
    assert data["summary"]["unchanged_pair_count"] == 15
    for category in ("escalated", "deescalated", "introduced", "resolved"):
        assert data[category] == []


# ---------- 错误码 ----------

def test_compare_missing_batch_404(client, batch_ids):
    base_id, target_id = batch_ids
    assert _compare(client, "no-such-batch", target_id).status_code == 404
    assert _compare(client, base_id, "no-such-batch").status_code == 404
    assert _compare(client, "no-such-batch", "neither").status_code == 404


# ---------- 配对键与分类逻辑单测 ----------

def test_pair_key_self_dimer_keeps_single_name():
    assert pair_key("F", "F", True) == ("F", "F")
    # 防御：同名即使未标记自二聚体也归并为单名键
    assert pair_key("F", "F", False) == ("F", "F")


def test_pair_key_heterodimer_ignores_order():
    assert pair_key("B", "A", False) == ("A", "B")
    assert pair_key("A", "B", False) == ("A", "B")


def test_compare_pair_risks_categories_and_sorting():
    base_pairs = [
        ("S", "S", True, "low"),
        ("B", "A", False, "low"),      # 键 (A,B)，提交顺序反置
        ("A", "C", False, "high"),
        ("A", "D", False, "medium"),
        ("A", "E", False, "low"),
    ]
    target_pairs = [
        ("S", "S", True, "high"),      # 自二聚体升级
        ("A", "B", False, "medium"),   # 升级
        ("C", "A", False, "low"),      # 降级（提交顺序反置）
        ("A", "D", False, "medium"),   # 不变
        ("A", "F", False, "low"),      # 新增
    ]
    result = compare_pair_risks(
        base_pairs=base_pairs,
        target_pairs=target_pairs,
        base_sequences={"A": "AAAA", "B": "CCCC", "C": "GGGG", "D": "TTTT", "E": "ACAC"},
        target_sequences={"A": "AAAA", "B": "CCCG", "C": "GGGG", "D": "TTTT", "F": "ACAC"},
    )

    assert [(d.name_a, d.name_b) for d in result.escalated] == [("A", "B"), ("S", "S")]
    assert [(d.name_a, d.name_b) for d in result.deescalated] == [("A", "C")]
    assert [(d.name_a, d.name_b) for d in result.introduced] == [("A", "F")]
    assert [(d.name_a, d.name_b) for d in result.resolved] == [("A", "E")]
    assert result.common_pair_count == 4
    assert result.unchanged_pair_count == 1

    # 等级与序列变化标记
    ab = result.escalated[0]
    assert (ab.base_risk_level, ab.target_risk_level) == ("low", "medium")
    assert ab.sequence_changed is True            # B 序列变化
    s_self = result.escalated[1]
    assert s_self.is_self_dimer and s_self.pair_key == "S"
    assert s_self.sequence_changed is False       # S 不在序列表中视为未变
    ac = result.deescalated[0]
    assert (ac.base_risk_level, ac.target_risk_level) == ("high", "low")
    assert ac.sequence_changed is False
    af = result.introduced[0]
    assert af.base_risk_level is None and af.target_risk_level == "low"
    ae = result.resolved[0]
    assert ae.target_risk_level is None and ae.base_risk_level == "low"


def test_compare_pair_risks_pair_key_display():
    result = compare_pair_risks(
        base_pairs=[],
        target_pairs=[("B", "A", False, "low"), ("S", "S", True, "high")],
        base_sequences={},
        target_sequences={"A": "AAAA", "B": "CCCC", "S": "GGGG"},
    )
    assert [d.pair_key for d in result.introduced] == ["A|B", "S"]
