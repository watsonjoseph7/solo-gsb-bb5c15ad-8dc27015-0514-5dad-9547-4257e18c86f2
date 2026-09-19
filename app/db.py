"""SQLite 持久化层：批次、引物、配对结果。"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone

from . import config

_SCHEMA = """
CREATE TABLE IF NOT EXISTS batches (
    id                       TEXT PRIMARY KEY,
    created_at               TEXT NOT NULL,
    continuous_threshold     INTEGER NOT NULL,
    three_prime_threshold    INTEGER NOT NULL,
    primer_count             INTEGER NOT NULL,
    pair_count               INTEGER NOT NULL,
    high_risk_count          INTEGER NOT NULL,
    medium_risk_count        INTEGER NOT NULL,
    low_risk_count           INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS primers (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id    TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    position    INTEGER NOT NULL,
    name        TEXT NOT NULL,
    sequence    TEXT NOT NULL,
    UNIQUE(batch_id, position)
);

CREATE TABLE IF NOT EXISTS pair_results (
    id                          INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id                    TEXT NOT NULL REFERENCES batches(id) ON DELETE CASCADE,
    pair_index                  INTEGER NOT NULL,
    name_top                    TEXT NOT NULL,
    name_bottom                 TEXT NOT NULL,
    is_self_dimer               INTEGER NOT NULL,
    alignments_checked          INTEGER NOT NULL,
    best_shift                  INTEGER NOT NULL,
    overlap_len                 INTEGER NOT NULL,
    max_complementary_run       INTEGER NOT NULL,
    top3_run                    INTEGER NOT NULL,
    bottom3_run                 INTEGER NOT NULL,
    terminal_run                INTEGER NOT NULL,
    involved_three_prime_ends   TEXT NOT NULL,
    risk_level                  TEXT NOT NULL,
    risk_reasons                TEXT NOT NULL,
    top_line                    TEXT NOT NULL,
    match_line                  TEXT NOT NULL,
    bottom_line                 TEXT NOT NULL,
    UNIQUE(batch_id, pair_index)
);

CREATE INDEX IF NOT EXISTS idx_pair_results_batch ON pair_results(batch_id);
CREATE INDEX IF NOT EXISTS idx_primers_batch ON primers(batch_id);
"""


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def init_db(reset: bool = False) -> None:
    """创建数据库目录与表结构。容器启动时执行。"""
    config.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(config.DB_PATH) as conn:
        conn.execute("PRAGMA journal_mode=WAL")
        if reset:
            conn.executescript(
                "DROP TABLE IF EXISTS pair_results;"
                "DROP TABLE IF EXISTS primers;"
                "DROP TABLE IF EXISTS batches;"
            )
        conn.executescript(_SCHEMA)


@contextmanager
def get_conn() -> Iterator[sqlite3.Connection]:
    conn = sqlite3.connect(config.DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys=ON")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def insert_batch(
    conn: sqlite3.Connection,
    *,
    batch_id: str,
    created_at: str,
    continuous_threshold: int,
    three_prime_threshold: int,
    primer_count: int,
    pair_count: int,
    risk_counts: dict[str, int],
    primers: list[tuple[int, str, str]],
    pair_rows: list[tuple],
) -> None:
    conn.execute(
        """INSERT INTO batches (id, created_at, continuous_threshold,
               three_prime_threshold, primer_count, pair_count,
               high_risk_count, medium_risk_count, low_risk_count)
           VALUES (?,?,?,?,?,?,?,?,?)""",
        (
            batch_id,
            created_at,
            continuous_threshold,
            three_prime_threshold,
            primer_count,
            pair_count,
            risk_counts.get("high", 0),
            risk_counts.get("medium", 0),
            risk_counts.get("low", 0),
        ),
    )
    conn.executemany(
        "INSERT INTO primers (batch_id, position, name, sequence) VALUES (?,?,?,?)",
        [(batch_id, pos, name, seq) for pos, name, seq in primers],
    )
    conn.executemany(
        """INSERT INTO pair_results (batch_id, pair_index, name_top, name_bottom,
               is_self_dimer, alignments_checked, best_shift, overlap_len,
               max_complementary_run, top3_run, bottom3_run, terminal_run,
               involved_three_prime_ends, risk_level, risk_reasons,
               top_line, match_line, bottom_line)
           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        [(batch_id, *row) for row in pair_rows],
    )


def get_batch(conn: sqlite3.Connection, batch_id: str) -> sqlite3.Row | None:
    return conn.execute(
        "SELECT * FROM batches WHERE id = ?", (batch_id,)
    ).fetchone()


def get_primers(conn: sqlite3.Connection, batch_id: str) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT position, name, sequence FROM primers "
        "WHERE batch_id = ? ORDER BY position",
        (batch_id,),
    ).fetchall()


def get_pairs(
    conn: sqlite3.Connection, batch_id: str, risky_only: bool = False
) -> list[sqlite3.Row]:
    sql = (
        "SELECT * FROM pair_results WHERE batch_id = ? "
        + ("AND risk_level != 'low' " if risky_only else "")
        + "ORDER BY pair_index"
    )
    return conn.execute(sql, (batch_id,)).fetchall()


def get_pair_metrics(
    conn: sqlite3.Connection, batch_id: str
) -> list[tuple[int, int]]:
    """只取阈值地图所需的两个标量指标：(最长连续互补, 3′端连续互补)。"""
    rows = conn.execute(
        "SELECT max_complementary_run, terminal_run FROM pair_results "
        "WHERE batch_id = ? ORDER BY pair_index",
        (batch_id,),
    ).fetchall()
    return [(r["max_complementary_run"], r["terminal_run"]) for r in rows]


def list_batches(
    conn: sqlite3.Connection, limit: int, offset: int
) -> list[sqlite3.Row]:
    return conn.execute(
        "SELECT * FROM batches ORDER BY created_at DESC, id "
        "LIMIT ? OFFSET ?",
        (limit, offset),
    ).fetchall()


def count_batches(conn: sqlite3.Connection) -> int:
    return conn.execute("SELECT COUNT(*) AS n FROM batches").fetchone()["n"]


def serialize_pair(row: sqlite3.Row) -> dict:
    return {
        "pair_index": row["pair_index"],
        "name_top": row["name_top"],
        "name_bottom": row["name_bottom"],
        "is_self_dimer": bool(row["is_self_dimer"]),
        "alignments_checked": row["alignments_checked"],
        "best_shift": row["best_shift"],
        "overlap_len": row["overlap_len"],
        "continuous_complementary_length": row["max_complementary_run"],
        "three_prime_run_length": row["terminal_run"],
        "top3_run_length": row["top3_run"],
        "bottom3_run_length": row["bottom3_run"],
        "involved_three_prime_ends": json.loads(row["involved_three_prime_ends"]),
        "risk_level": row["risk_level"],
        "risk_reasons": json.loads(row["risk_reasons"]),
        "alignment": {
            "top": row["top_line"],
            "match": row["match_line"],
            "bottom": row["bottom_line"],
        },
    }
