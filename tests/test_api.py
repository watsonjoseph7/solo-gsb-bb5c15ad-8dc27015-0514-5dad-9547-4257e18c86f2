"""API 端到端测试（FastAPI TestClient + 临时 SQLite）。"""

from __future__ import annotations

import os

import pytest


@pytest.fixture(scope="session")
def client(tmp_path_factory):
    # 必须在导入 app.config 之前指定临时数据库路径
    tmp_dir = tmp_path_factory.mktemp("db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},
        {"name": "R", "sequence": "GGGGGTTTTT"},
        {"name": "degenerate", "sequence": "ACGTNRYSWKMBDHV"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


def test_health(client):
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_create_batch_pair_count_and_fields(client):
    resp = client.post("/api/batches", json=PAYLOAD)
    assert resp.status_code == 201, resp.text
    data = resp.json()

    # 3 自二聚体 + 3 异二聚体
    assert data["primer_count"] == 3
    assert data["pair_count"] == 6
    assert len(data["pairs"]) == 6

    self_pairs = [p for p in data["pairs"] if p["is_self_dimer"]]
    hetero_pairs = [p for p in data["pairs"] if not p["is_self_dimer"]]
    assert len(self_pairs) == 3
    assert len(hetero_pairs) == 3

    for p in data["pairs"]:
        assert p["alignments_checked"] >= 1
        assert p["risk_level"] in {"high", "medium", "low"}
        assert p["risk_reasons"]
        assert len(p["alignment"]["top"]) == len(p["alignment"]["bottom"])
        assert "continuous_complementary_length" in p
        assert "three_prime_run_length" in p
        assert "best_shift" in p

    # F/R 是典型 3′端二聚体
    fr = next(p for p in hetero_pairs if p["name_top"] == "F" and p["name_bottom"] == "R")
    assert fr["risk_level"] == "high"
    assert fr["continuous_complementary_length"] == 5
    assert fr["three_prime_run_length"] == 5
    assert fr["top3_run_length"] == 5
    assert fr["bottom3_run_length"] == 5
    assert set(fr["involved_three_prime_ends"]) == {"F:3′", "R:3′"}


def test_get_batch_persistence_and_risky_filter(client):
    created = client.post("/api/batches", json=PAYLOAD).json()
    bid = created["id"]

    full = client.get(f"/api/batches/{bid}")
    assert full.status_code == 200
    assert full.json()["pair_count"] == 6

    risky = client.get(f"/api/batches/{bid}?risky_only=true").json()
    assert all(p["risk_level"] in {"high", "medium"} for p in risky["pairs"])
    assert len(risky["pairs"]) <= 6

    missing = client.get("/api/batches/does-not-exist")
    assert missing.status_code == 404


def test_list_batches(client):
    listing = client.get("/api/batches").json()
    assert listing["total"] >= 2
    assert len(listing["batches"]) == listing["total"]


def test_validation_errors(client):
    # 非法 IUPAC
    bad = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "x", "sequence": "ACGX"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 2,
        },
    )
    assert bad.status_code == 422

    # 阈值必须 >=1
    bad_thr = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "x", "sequence": "ACGT"}],
            "continuous_threshold": 0,
            "three_prime_threshold": 2,
        },
    )
    assert bad_thr.status_code == 422

    # 空批次
    empty = client.post(
        "/api/batches",
        json={
            "primers": [],
            "continuous_threshold": 3,
            "three_prime_threshold": 2,
        },
    )
    assert empty.status_code == 422

    # 名称重复
    dup = client.post(
        "/api/batches",
        json={
            "primers": [
                {"name": "x", "sequence": "ACGT"},
                {"name": "x", "sequence": "TTTT"},
            ],
            "continuous_threshold": 3,
            "three_prime_threshold": 2,
        },
    )
    assert dup.status_code == 422
    assert "重复" in dup.json()["detail"]


def test_iupac_degenerate_complementarity(client):
    payload = {
        "primers": [
            {"name": "a", "sequence": "RRRR"},   # A/G
            {"name": "b", "sequence": "YYYY"},   # C/T，与 R 互补
        ],
        "continuous_threshold": 4,
        "three_prime_threshold": 3,
    }
    data = client.post("/api/batches", json=payload).json()
    hetero = next(
        p for p in data["pairs"]
        if not p["is_self_dimer"] and p["name_top"] == "a"
    )
    assert hetero["continuous_complementary_length"] == 4
    assert hetero["risk_level"] == "high"


def test_single_primer_batch_is_self_dimer_only(client):
    data = client.post(
        "/api/batches",
        json={
            "primers": [{"name": "solo", "sequence": "AAAATTTT"}],
            "continuous_threshold": 3,
            "three_prime_threshold": 3,
        },
    ).json()
    assert data["primer_count"] == 1
    assert data["pair_count"] == 1
    assert data["pairs"][0]["is_self_dimer"] is True
    # 取回的批次与提交响应一致
    fetched = client.get(f"/api/batches/{data['id']}").json()
    assert fetched["pairs"][0]["name_top"] == "solo"
