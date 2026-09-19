"""简单命令行：initdb 初始化/校验数据库。"""

from __future__ import annotations

import sys

from . import config, db


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "initdb":
        db.init_db(reset=config.RESET_ON_START)
        print(f"database ready at {config.DB_PATH}")
        return 0
    print("usage: python -m app.cli initdb", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
