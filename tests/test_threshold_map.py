"""阈值地图接口端到端测试与二维差分聚合单测。"""

from __future__ import annotations

import os
import random

import pytest

from app.schemas import ThresholdMapRequest
from app.threshold_map import build_threshold_map


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("threshold_map_db")
    os.environ["PRIMER_DB_PATH"] = str(tmp_dir / "test.db")

    from app import db, main

    db.init_db()
    from fastapi.testclient import TestClient

    with TestClient(main.app) as c:
        yield c


# 6 引物批次（21 个配对），指标已经 dimer 算法核实：
#   (0,0)×6  (4,0)×2  (2,2)×1  (3,3)×3  (4,4)×4  (5,5)×2  (8,8)×2  (13,13)×1
PAYLOAD = {
    "primers": [
        {"name": "F", "sequence": "GGGGGAAAAA"},
        {"name": "R", "sequence": "GGGGGTTTTT"},
        {"name": "X", "sequence": "AAAAAAAAAA"},
        {"name": "H", "sequence": "AAAATTTT"},
        {"name": "D", "sequence": "ACGTNRYSWKMBDHV"},
        {"name": "S", "sequence": "CCCCGGGG"},
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3,
}


@pytest.fixture(scope="module")
def batch_id(client):
    return client.post("/api/batches", json=PAYLOAD).json()["id"]


def _map(client, batch_id, **over):
    body = {
        "continuous_start": 1,
        "continuous_end": 8,
        "three_prime_start": 1,
        "three_prime_end": 8,
    }
    body.update(over)
    return client.post(f"/api/batches/{batch_id}/threshold-map", json=body)


def _brute_force(metrics, ct, tt):
    """现有判定规则的逐配对朴素实现：high 优先、其次 medium、否则 low。"""
    high = sum(1 for _c, t in metrics if t >= tt)
    medium = sum(1 for c, t in metrics if t < tt and c >= ct)
    return high, medium, len(metrics) - high - medium


def _saved_metrics(client, batch_id):
    pairs = client.get(f"/api/batches/{batch_id}").json()["pairs"]
    return [
        (p["continuous_complementary_length"], p["three_prime_run_length"])
        for p in pairs
    ]


# ---------- 聚合正确性 ----------

def test_every_cell_matches_brute_force(client, batch_id):
    metrics = _saved_metrics(client, batch_id)
    assert len(metrics) == 21

    resp = _map(client, batch_id)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["batch_id"] == batch_id
    assert data["pair_count"] == 21
    assert data["cell_count"] == 64
    assert len(data["cells"]) == 64

    for cell in data["cells"]:
        ct, tt = cell["continuous_threshold"], cell["three_prime_threshold"]
        high, medium, low = _brute_force(metrics, ct, tt)
        assert (cell["high_count"], cell["medium_count"], cell["low_count"]) == (
            high,
            medium,
            low,
        ), f"cell ({ct}, {tt})"
        # 每个单元三级计数之和恒等于配对总数
        assert high + medium + low == 21


def test_known_spot_values(client, batch_id):
    cells = {
        (c["continuous_threshold"], c["three_prime_threshold"]): c
        for c in _map(client, batch_id, three_prime_end=10).json()["cells"]
    }
    # tt=1 时 3′端指标 >=1 的配对全部 high（(0,0)×6 与 (4,0)×2 共 8 对例外）
    assert cells[(1, 1)]["high_count"] == 13
    # tt=9：仅 (13,13) 一对 high；(8,8)×2 落入 medium（ct<=8）
    assert cells[(8, 9)]["high_count"] == 1
    assert cells[(8, 9)]["medium_count"] == 2
    # ct=6、tt=6：(5,5) 及以下不再 high/medium；high 为 (8,8)×2+(13,13)
    assert cells[(6, 6)]["high_count"] == 3
    assert cells[(6, 6)]["medium_count"] == 0
    # (4,0)×2 对：ct<=4 且 tt>=1 时永远 medium（t=0 不可能 high）
    assert cells[(4, 1)]["medium_count"] == 2
    assert cells[(5, 1)]["medium_count"] == 0


def test_batch_threshold_cell_marked_and_matches_batch_counts(client, batch_id):
    data = _map(client, batch_id).json()

    marked = [c for c in data["cells"] if c["is_batch_threshold"]]
    assert len(marked) == 1
    cell = marked[0]
    assert cell["continuous_threshold"] == 4
    assert cell["three_prime_threshold"] == 3

    # 与批次落库的三级计数一致（同一判定规则、同一批指标）
    batch = client.get(f"/api/batches/{batch_id}").json()
    assert cell["high_count"] == batch["high_risk_count"] == 12
    assert cell["medium_count"] == batch["medium_risk_count"] == 2
    assert cell["low_count"] == batch["low_risk_count"] == 7

    assert data["batch_continuous_threshold"] == 4
    assert data["batch_three_prime_threshold"] == 3


def test_batch_threshold_outside_range_marks_nothing(client, batch_id):
    # 批次原阈值 (4, 3) 不在请求范围内时不应有任何标记单元
    data = _map(
        client, batch_id,
        continuous_start=5, continuous_end=8,
        three_prime_start=4, three_prime_end=8,
    ).json()
    assert data["cell_count"] == 4 * 5
    assert all(not c["is_batch_threshold"] for c in data["cells"])


def test_cells_sorted_by_three_prime_then_continuous(client, batch_id):
    cells = _map(client, batch_id).json()["cells"]
    keys = [(c["three_prime_threshold"], c["continuous_threshold"]) for c in cells]
    assert keys == sorted(keys)
    assert keys[0] == (1, 1)
    assert keys[-1] == (8, 8)
    # 每个阈值组合恰好出现一次
    assert len(set(keys)) == len(keys)


def test_single_cell_grid(client, batch_id):
    data = _map(
        client, batch_id,
        continuous_start=4, continuous_end=4,
        three_prime_start=3, three_prime_end=3,
    ).json()
    assert data["cell_count"] == 1
    cell = data["cells"][0]
    assert cell["is_batch_threshold"] is True
    assert (cell["high_count"], cell["medium_count"], cell["low_count"]) == (
        12, 2, 7,
    )


def test_full_range_grid_accepted(client, batch_id):
    # 200×200 = 40000 恰好达到组合数上限，必须接受
    resp = _map(
        client, batch_id,
        continuous_start=1, continuous_end=200,
        three_prime_start=1, three_prime_end=200,
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["cell_count"] == 40000
    assert len(data["cells"]) == 40000
    # 阈值超过最大指标（13）后全部配对 low
    corner = next(
        c for c in data["cells"]
        if c["continuous_threshold"] == 14 and c["three_prime_threshold"] == 14
    )
    assert (corner["high_count"], corner["medium_count"], corner["low_count"]) == (
        0, 0, 21,
    )


# ---------- 只读保证 ----------

def test_no_alignment_recompute(client, batch_id, monkeypatch):
    import app.services as services

    calls = []
    original = services.analyze_pair

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(services, "analyze_pair", spy)
    resp = _map(client, batch_id)
    assert resp.status_code == 200
    assert calls == []


def test_batch_not_modified(client, batch_id):
    before = client.get(f"/api/batches/{batch_id}").json()
    listing_before = client.get("/api/batches").json()["total"]

    resp = _map(client, batch_id)
    assert resp.status_code == 200

    after = client.get(f"/api/batches/{batch_id}").json()
    assert after == before
    assert client.get("/api/batches").json()["total"] == listing_before


# ---------- 错误码 ----------

def test_batch_not_found_404(client):
    resp = _map(client, "no-such-batch")
    assert resp.status_code == 404
    assert "不存在" in resp.json()["detail"]


def test_inverted_ranges_422(client, batch_id):
    resp = _map(client, batch_id, continuous_start=5, continuous_end=2)
    assert resp.status_code == 422
    resp = _map(client, batch_id, three_prime_start=8, three_prime_end=3)
    assert resp.status_code == 422


def test_threshold_bounds_422(client, batch_id):
    assert _map(client, batch_id, continuous_start=0).status_code == 422
    assert _map(client, batch_id, three_prime_end=0).status_code == 422
    assert _map(client, batch_id, continuous_end=201).status_code == 422
    assert _map(client, batch_id, three_prime_start=201).status_code == 422


def test_combination_limit_validator():
    # 40000 恰好允许；超出（绕过字段上界构造）必须拒绝
    ok = ThresholdMapRequest(
        continuous_start=1, continuous_end=200,
        three_prime_start=1, three_prime_end=200,
    )
    assert ok._check_ranges() is ok

    too_many = ThresholdMapRequest.model_construct(
        continuous_start=1, continuous_end=200,
        three_prime_start=1, three_prime_end=201,
    )
    with pytest.raises(ValueError, match="40000"):
        too_many._check_ranges()


# ---------- 二维差分聚合单测 ----------

def test_build_threshold_map_matches_brute_force_random():
    rng = random.Random(20260919)
    for _ in range(50):
        metrics = [(rng.randint(0, 30), rng.randint(0, 30)) for _ in range(80)]
        c_start = rng.randint(1, 10)
        c_end = rng.randint(c_start, 40)
        t_start = rng.randint(1, 10)
        t_end = rng.randint(t_start, 40)
        grid = build_threshold_map(
            metrics,
            continuous_start=c_start, continuous_end=c_end,
            three_prime_start=t_start, three_prime_end=t_end,
        )
        assert len(grid) == t_end - t_start + 1
        assert all(len(row) == c_end - c_start + 1 for row in grid)
        for r, tt in enumerate(range(t_start, t_end + 1)):
            for c, ct in enumerate(range(c_start, c_end + 1)):
                assert grid[r][c] == _brute_force(metrics, ct, tt)


def test_build_threshold_map_edge_cases():
    # 空指标：全部单元为 (0, 0, 0)
    grid = build_threshold_map(
        [], continuous_start=1, continuous_end=3,
        three_prime_start=1, three_prime_end=3,
    )
    assert grid == [[(0, 0, 0)] * 3 for _ in range(3)]

    # 指标恰好等于阈值边界：>= 语义
    grid = build_threshold_map(
        [(5, 5)], continuous_start=5, continuous_end=6,
        three_prime_start=5, three_prime_end=6,
    )
    # tt=5：high（t>=5）；tt=6：不再 high，ct=5 时 medium，ct=6 时 low
    assert grid[0] == [(1, 0, 0), (1, 0, 0)]
    assert grid[1] == [(0, 1, 0), (0, 0, 1)]

    # 指标超出网格范围：按边界截断后仍正确
    grid = build_threshold_map(
        [(200, 200), (0, 0)], continuous_start=3, continuous_end=4,
        three_prime_start=3, three_prime_end=4,
    )
    assert grid == [[(1, 0, 1), (1, 0, 1)], [(1, 0, 1), (1, 0, 1)]]
