import os
import tempfile

# 在任何 app 模块被收集导入前，指向临时数据库
_TMP = tempfile.mkdtemp(prefix="primer_test_")
os.environ.setdefault("PRIMER_DB_PATH", os.path.join(_TMP, "test.db"))
