#!/bin/sh
# 容器启动入口：初始化 SQLite 数据库（建表）后启动 API
set -e

echo "[entrypoint] initializing SQLite database at ${PRIMER_DB_PATH:-/data/primer_dimer.db}"
python3 -m app.cli initdb

echo "[entrypoint] starting API on ${PRIMER_HOST:-0.0.0.0}:${PRIMER_PORT:-8000}"
exec python3 -m uvicorn app.main:app \
    --host "${PRIMER_HOST:-0.0.0.0}" \
    --port "${PRIMER_PORT:-8000}"
