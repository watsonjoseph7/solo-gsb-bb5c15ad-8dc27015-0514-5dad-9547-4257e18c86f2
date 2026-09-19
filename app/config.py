from __future__ import annotations

import os
from pathlib import Path


def _env_bool(name: str, default: str = "0") -> bool:
    return os.getenv(name, default).strip().lower() in {"1", "true", "yes", "on"}


# SQLite 数据库文件路径（目录不存在时自动创建）
DB_PATH: Path = Path(os.getenv("PRIMER_DB_PATH", "/data/primer_dimer.db"))

# 服务监听地址与端口
HOST: str = os.getenv("PRIMER_HOST", "0.0.0.0")
PORT: int = int(os.getenv("PRIMER_PORT", "8000"))

# 输入上限：避免畸形超大批次拖垮服务
MAX_PRIMERS_PER_BATCH: int = int(os.getenv("PRIMER_MAX_PRIMERS", "500"))
MAX_PRIMER_LENGTH: int = int(os.getenv("PRIMER_MAX_LENGTH", "200"))

# 数据库初始化时是否重建表结构（容器首次启动自动初始化，无需开启）
RESET_ON_START: bool = _env_bool("PRIMER_DB_RESET", "0")

# 阈值地图接口的固定约束（接口规格，非部署可调）：
# 阈值起止范围 1～200，单次请求的阈值组合数不超过 40000
THRESHOLD_MAP_MAX_THRESHOLD: int = 200
THRESHOLD_MAP_MAX_CELLS: int = 40000
