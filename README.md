# PCR 引物二聚体筛查 API

基于 **Python + FastAPI + SQLite** 的引物二聚体（primer dimer）批量筛查服务。
调用方提交一批带名称的 IUPAC DNA 引物及两个阈值（连续互补阈值、3′端连续互补阈值），
服务枚举每对引物**自二聚体**与**两两异二聚体的全部错位对齐**，按 IUPAC 简并碱基
**集合语义**判断碱基互补，返回每对的最危险对齐、最长连续互补长度、涉及任一 3′端的
连续长度及中文风险原因；同分时按固定偏移规则取唯一结果。批次与结果持久化到
SQLite，支持事后查询。对已保存批次还可通过 `POST /api/batches/{id}/pool-plan`
基于异二聚体结果构建冲突图，按 DSatur 风格确定性贪心把引物分配到若干引物池；
通过 `POST /api/batches/{id}/replacement-ranking` 在**内存**中模拟用候选序列
替换某条引物，按 high 配对数、最大 3′连续互补、最大连续互补三个指标做
非支配分层排名（不落库）；
通过 `POST /api/batches/{id}/augmentation-plan` 在**内存**中对 2～18 条候选
引物复用现有配对算法计算自二聚体、与原批次及候选两两异二聚体，在
high/medium 冲突图上用**确定性分支限界**精确选出不超过接纳上限的最大无冲突
候选集合（不落库）；
通过 `POST /api/batches/{id}/threshold-map` 只读取已保存配对的两个标量指标
（最长连续互补、3′端连续互补），用**二维差分**在请求的阈值网格上聚合出每个
阈值组合下的 high/medium/low 配对数量（不重算对齐、不修改批次）；
通过 `GET /api/batches/{base_id}/compare/{target_id}` 按引物名称构造稳定配对键，
对比两个批次的配对风险升降（escalated / deescalated / introduced / resolved）。

---

## 一条命令启动

```bash
docker compose up --build -d
```

启动后访问：

| 用途 | 地址 |
| --- | --- |
| 交互式 API 文档（Swagger UI） | http://localhost:8000/docs |
| ReDoc 文档 | http://localhost:8000/redoc |
| 健康检查 | http://localhost:8000/health |
| 批次接口 | http://localhost:8000/api/batches |
| 引物分池 | `POST http://localhost:8000/api/batches/{id}/pool-plan` |
| 替换候选排名 | `POST http://localhost:8000/api/batches/{id}/replacement-ranking` |
| 候选增量筛选 | `POST http://localhost:8000/api/batches/{id}/augmentation-plan` |
| 阈值扫描地图 | `POST http://localhost:8000/api/batches/{id}/threshold-map` |
| 批次风险对比 | `GET http://localhost:8000/api/batches/{base_id}/compare/{target_id}` |

> 监听端口 **8000**。SQLite 数据库文件保存在命名卷 `primer-data` 内
> （容器路径 `/data/primer_dimer.db`）。镜像构建期安装依赖、容器启动期自动建库
> （`docker-entrypoint.sh` → `python -m app.cli initdb` → `uvicorn`）。

停止与查看日志：

```bash
docker compose logs -f
docker compose down            # 保留数据卷
docker compose down -v         # 同时删除数据库卷
```

---

## 必要配置（环境变量，均有默认值，不配也能直接启动）

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PRIMER_DB_PATH` | `/data/primer_dimer.db` | SQLite 数据库文件路径 |
| `PRIMER_HOST` | `0.0.0.0` | 监听地址 |
| `PRIMER_PORT` | `8000` | 监听端口（容器内；改宿主机端口请编辑 `docker-compose.yml` 的 `ports`） |
| `PRIMER_MAX_PRIMERS` | `500` | 单批次引物数量上限 |
| `PRIMER_MAX_LENGTH` | `200` | 单条引物长度上限 |
| `PRIMER_DB_RESET` | `0` | 设为 `1` 时启动清空并重建表（仅用于重置数据） |

---

## API 概览

### `POST /api/batches` — 提交批次（同步计算后返回完整结果）

```bash
curl -s -X POST http://localhost:8000/api/batches \
  -H 'Content-Type: application/json' \
  -d '{
    "primers": [
      {"name": "F", "sequence": "GGGGGAAAAA"},
      {"name": "R", "sequence": "GGGGGTTTTT"}
    ],
    "continuous_threshold": 4,
    "three_prime_threshold": 3
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `primers[].name` | 引物名称，非空，批次内唯一 |
| `primers[].sequence` | 5′→3′ IUPAC DNA 序列，接受 `A T C G R Y S W K M B D H V N`（大小写、首尾空白自动规范化；含 `-`、空格等非 IUPAC 字符返回 422） |
| `continuous_threshold` | 连续互补长度阈值（≥1），达到即提示双链结合稳定 |
| `three_prime_threshold` | 3′端连续互补长度阈值（≥1），达到即提示聚合酶延伸风险，风险等级 `high` |

返回（节选）：

```json
{
  "id": "f98cdd27-...",
  "primer_count": 2,
  "pair_count": 3,
  "high_risk_count": 1, "medium_risk_count": 0, "low_risk_count": 2,
  "pairs": [
    {
      "pair_index": 2,
      "name_top": "F", "name_bottom": "R", "is_self_dimer": false,
      "alignments_checked": 19,
      "best_shift": 5,
      "overlap_len": 10,
      "continuous_complementary_length": 5,
      "three_prime_run_length": 5,
      "top3_run_length": 5, "bottom3_run_length": 5,
      "involved_three_prime_ends": ["F:3′", "R:3′"],
      "risk_level": "high",
      "risk_reasons": [
        "3′端连续互补达到阈值，存在聚合酶延伸风险：top 3′端连续互补 5 bp，bottom 3′端连续互补 5 bp",
        "最长连续互补 5 bp 达到阈值 4 bp，双链结合稳定"
      ],
      "alignment": {
        "top":    "GGGGGAAAAA     ",
        "match":  "     |||||     ",
        "bottom": "     TTTTTGGGGG"
      }
    }
  ]
}
```

### `GET /api/batches/{id}` — 批次结果查询

* `GET /api/batches/{id}` 返回该批次全部配对；
* `GET /api/batches/{id}?risky_only=true` 只返回 `high`/`medium` 风险配对；
* 不存在的批次返回 `404`。

### `GET /api/batches?limit=50&offset=0` — 批次列表（按创建时间倒序）

### `POST /api/batches/{id}/pool-plan` — 对已保存批次做引物分池

对**已保存批次**的异二聚体筛查结果构建冲突图并分池：

* **自二聚体忽略**，只取异二聚体配对；风险等级命中 `isolated_risk_levels`
  （`high` / `medium` / `low`，可多选、自动去重）的异二聚体两端引物之间连一条冲突边；
* 贪心选择下一个引物，依次比较 **① 已见邻池数**（相邻引物已占据的不同池数，
  DSatur saturation）、**② 冲突度**（冲突图中的边数）、**③ 批次位置**（靠前者优先）；
* 贪心选择放置池，依次比较 **① 新增冲突数**（放入后新产生的同池冲突边数，少者优先）、
  **② 池内人数**（少者优先）、**③ 池号**（小者优先）；
* 容量允许（`pool_count × pool_capacity ≥ 引物数`）时**全部引物必被分配**；
  容量总和不足返回 `422`；批次不存在返回 `404`。全部比较均为全序，相同输入结果稳定可复现。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../pool-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "pool_count": 3,
    "pool_capacity": 8,
    "isolated_risk_levels": ["high", "medium"]
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `pool_count` | 池数量，≥1 |
| `pool_capacity` | 每池容量（每池最多引物数），≥1 |
| `isolated_risk_levels` | 需隔离的异二聚体风险等级列表，取值 `high`/`medium`/`low`，至少 1 个；重复值去重 |

返回示例：

```json
{
  "batch_id": "f98cdd27-...",
  "pool_count": 3,
  "pool_capacity": 8,
  "isolated_risk_levels": ["high", "medium"],
  "primer_count": 3,
  "assigned_count": 3,
  "total_conflicts": 0,
  "pools": [
    {
      "pool": 0,
      "capacity": 8,
      "members": ["R"],
      "member_positions": [1],
      "conflict_count": 0
    }
  ],
  "conflicts": [
    {
      "pool": 0,
      "position_a": 1, "position_b": 4,
      "name_a": "R", "name_b": "P5",
      "risk_level": "high"
    }
  ]
}
```

* `pools[].members` / `member_positions` 为各池成员，按批次位置升序；
  `conflict_count` 为该池**内部**同池冲突边数；
* `conflicts` 平铺全部池的同池冲突边（按池号、两端批次位置升序），
  长度恒等于 `total_conflicts`，也等于各池 `conflict_count` 之和；
* 容量允许但冲突图无法完全隔离时（如池数少于色数），算法仍分配全部引物，
  不可避免的冲突通过 `conflicts` 如实返回。

### `POST /api/batches/{id}/replacement-ranking` — 候选替换引物模拟与非支配分层

对**已保存批次**中的一条目标引物，提交 **2～50 条带唯一标签的候选 IUPAC 序列**，
服务在**内存**中逐条模拟“用候选替换目标引物”后的二聚体情况并排名：

* **只重算受影响配对**：每个候选仅复用现有 `analyze_pair` 二聚体算法计算
  ① 候选的**自二聚体**、② 候选与**其余每条引物**的异二聚体；目标以外引物
  彼此之间的既有配对**不重算**；
* **沿用批次原阈值**：风险等级（`high`/`medium`/`low`）按批次提交时的
  `continuous_threshold`、`three_prime_threshold` 判定，无需在请求中重复给出；
* **三个汇总指标（均越小越优）**：
  1. `high_pair_count` —— 受影响配对中风险等级为 `high` 的数量（含候选自二聚体）；
  2. `max_three_prime_run` —— 受影响配对中的最大 3′连续互补长度；
  3. `max_continuous_run` —— 受影响配对中的最长连续互补长度；
* **非支配分层**：候选 A 在三个指标上均不劣于 B 且至少一维严格更优，则 A
  支配 B；全体非支配候选为第 1 层，移除后再取非支配集合为第 2 层，依此类推。
  **同层候选按请求输入顺序稳定排序**（指标完全相同的候选必然同层）；
* **不落库**：纯内存模拟，不修改批次、引物与已保存配对，也不会新建批次；
* 批次不存在返回 `404`；目标引物名不存在、候选标签重复（去除首尾空白后比较）、
  候选数量超出 2～50 或序列含非法 IUPAC 字符返回 `422`。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../replacement-ranking \
  -H 'Content-Type: application/json' \
  -d '{
    "target_name": "F",
    "candidates": [
      {"label": "keep", "sequence": "GGGGGAAAAA"},
      {"label": "better-3prime", "sequence": "GGGGGAAAA"},
      {"label": "hairpin-risk", "sequence": "CCCCCGGGGG"}
    ]
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `target_name` | 待替换引物名称，必须为批次内已存在的引物名（不存在返回 422） |
| `candidates` | 2～50 条候选；`label` 非空、请求内唯一（去首尾空白后比较）；`sequence` 为 5′→3′ IUPAC DNA 序列，校验规则与批次提交一致 |

返回示例：

```json
{
  "batch_id": "f98cdd27-...",
  "target_name": "F",
  "continuous_threshold": 4,
  "three_prime_threshold": 3,
  "primer_count": 2,
  "candidate_count": 3,
  "candidates": [
    {
      "label": "better-3prime",
      "sequence": "GGGGGAAAA",
      "layer": 1,
      "rank": 1,
      "high_pair_count": 1,
      "max_three_prime_run": 4,
      "max_continuous_run": 4,
      "high_pairs": [ { "pair_index": 2, "name_top": "F", "name_bottom": "R", "...": "与批次配对结果同构" } ],
      "max_three_prime_pairs": [ { "pair_index": 2, "...": "..." } ],
      "max_continuous_pairs": [ { "pair_index": 2, "...": "..." } ]
    },
    {
      "label": "keep",
      "sequence": "GGGGGAAAAA",
      "layer": 2, "rank": 2,
      "high_pair_count": 1,
      "max_three_prime_run": 5,
      "max_continuous_run": 5,
      "high_pairs": [ { "pair_index": 2, "...": "..." } ],
      "max_three_prime_pairs": [ { "pair_index": 2, "...": "..." } ],
      "max_continuous_pairs": [ { "pair_index": 2, "...": "..." } ]
    },
    {
      "label": "hairpin-risk",
      "sequence": "CCCCCGGGGG",
      "layer": 3, "rank": 3,
      "high_pair_count": 1,
      "max_three_prime_run": 10,
      "max_continuous_run": 10,
      "high_pairs": [ { "pair_index": 0, "name_top": "F", "name_bottom": "F", "is_self_dimer": true, "...": "..." } ],
      "max_three_prime_pairs": [ { "pair_index": 0, "...": "..." } ],
      "max_continuous_pairs": [ { "pair_index": 0, "...": "..." } ]
    }
  ]
}
```

* `layer` 从 1 起（第 1 层即 Pareto 最优集合），`rank` 从 1 起、按层递增；
  候选的最终顺序为先层级、再输入顺序，稳定可复现；
* `high_pairs` / `max_three_prime_pairs` / `max_continuous_pairs` 中的每个配对
  与 `GET /api/batches/{id}` 的配对结果字段**完全同构**（含 `pair_index`、
  `risk_level`、对齐可视化等）；配对中的引物名沿用目标引物名（候选不入库，
  不产生新名称）；
* `pair_index` 与原批次落库位置一致（自二聚体占据目标引物位置，异二聚体按
  提交顺序枚举），便于与原批次结果按位置对照；
* 达到最大值的配对**全部**列出（可能有并列多个，包括候选自二聚体）；
  最大 3′连续互补为 0 时，所有该指标为 0 的受影响配对都会列入
  `max_three_prime_pairs`。

### `POST /api/batches/{id}/augmentation-plan` — 候选引物增量筛选（冲突图精确求解，不落库）

向**已保存批次**追加候选引物时，提交 **2～18 条名称唯一的候选 IUPAC 序列**与
最多接纳数 `max_selected`（1～18），服务在**内存**中计算候选的二聚体风险并
精确给出可同时接纳的最大无冲突候选集合：

* **复用现有配对算法、沿用批次原阈值**：对每个候选复用 `analyze_pair` 计算
  ① 候选的**自二聚体**、② 候选与**每条原批次引物**的异二聚体、③ **候选两两**
  异二聚体；风险等级按批次提交时的 `continuous_threshold`、
  `three_prime_threshold` 判定，无需在请求中重复给出；
* **冲突图**：候选为顶点。候选自二聚体或候选与原批次引物的配对达到
  `high`/`medium` 时，该候选记**自环**（原批次引物固定在场，等价于候选被固定
  禁止）；候选两两的 `high`/`medium` 配对记普通无向边（`low` 配对不产生边）；
* **确定性分支限界精确解**：在冲突图上求**不超过 `max_selected`** 的最大独立集，
  按候选下标升序“包含优先 / 排除其后”枚举（以贪心团覆盖做独立数上界剪枝），
  故同规模时**按请求顺序字典序最小**的解为唯一返回，结果稳定可复现；
* **淘汰原因**：
  * `risk` —— 存在阻断配对，`blocking_pairs` 列出与**原批次引物**或**入选候选**
    的全部阻断配对（含候选自二聚体），每个配对与 `GET /api/batches/{id}` 的
    配对结果字段同构，并用 `counterpart_kind`（`self`/`batch`/`candidate`）与
    `counterpart_name` 指明配对另一端；
  * `capacity`（`limited_by_capacity: true`）—— 候选本身无任何阻断配对，仅因
    不超过 `max_selected` 的数量上限未入选，`blocking_pairs` 为空；
* **不落库**：纯内存计算，不修改批次、引物与已保存配对，也不新建批次；
* 批次不存在返回 `404`；候选数量超出 2～18、`max_selected` 超出 1～18、
  候选名称为空/重复（去首尾空白后比较）、**与批次内已有引物重名**、或序列含
  非法 IUPAC 字符返回 `422`。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../augmentation-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "max_selected": 10,
    "candidates": [
      {"name": "newF", "sequence": "ACGTACGTAA"},
      {"name": "newR", "sequence": "TTTTGGGGCC"}
    ]
  }'
```

返回示例：

```json
{
  "batch_id": "f98cdd27-...",
  "continuous_threshold": 4,
  "three_prime_threshold": 3,
  "max_selected": 10,
  "candidate_count": 2,
  "selected_count": 1,
  "rejected_count": 1,
  "selected": [
    {"position": 0, "name": "newF", "sequence": "ACGTACGTAA"}
  ],
  "rejected": [
    {
      "position": 1,
      "name": "newR",
      "sequence": "TTTTGGGGCC",
      "reason": "risk",
      "limited_by_capacity": false,
      "blocking_pairs": [
        {
          "pair_index": 3,
          "name_top": "F", "name_bottom": "newR", "is_self_dimer": false,
          "risk_level": "high",
          "continuous_complementary_length": 4,
          "three_prime_run_length": 4,
          "counterpart_kind": "batch",
          "counterpart_name": "F",
          "alignment": { "...": "与批次配对结果同构" }
        }
      ]
    }
  ]
}
```

* `selected` / `rejected` 均按请求顺序（`position` 即候选在请求中的序号，0 起）；
* `blocking_pairs[].pair_index` 为本次内存枚举内的稳定序号：`[0, m)` 为候选
  自二聚体，其后按候选、批次位置排列候选-原批次配对，最后按 `(i,j)` 行枚举
  候选两两配对；
* 被自身自二聚体淘汰时 `counterpart_kind: "self"`；与入选候选冲突时
  `counterpart_kind: "candidate"`（只会列出最终入选的候选；与其他落选候选的
  冲突不构成阻断）。

### `POST /api/batches/{id}/threshold-map` — 阈值扫描地图

对**已保存批次**做一次“阈值敏感性扫描”：请求给出连续互补阈值与 3′端阈值的
整数起止范围（均 **1～200**，闭区间），服务返回网格中**每个阈值组合**下判定为
`high` / `medium` / `low` 的配对数量，用于评估阈值取值对风险结论的影响：

* **只读已保存指标**：仅使用每个配对已落库的两个标量——最长连续互补长度
  `continuous_complementary_length` 与 3′端连续互补长度
  `three_prime_run_length`；**不重算任何对齐**，也**不修改批次**（纯读操作）；
* **判定规则与批次计算完全一致**：`3′端连续 ≥ 3′端阈值` → `high`（最优先）；
  否则 `最长连续互补 ≥ 连续阈值` → `medium`；否则 `low`；
* **二维差分聚合**：每个配对对 high、medium 计数各贡献一个轴对齐矩形阈值
  区域，用二维差分数组（imos）打角点标记后一次前缀和物化整图，复杂度
  O(配对数 + 单元数)，**不逐配对遍历网格**；
* **组合数上限**：`(continuous_end - continuous_start + 1) ×
  (three_prime_end - three_prime_start + 1) ≤ 40000`（恰好允许完整的
  200×200 网格），超限或起点大于终点返回 `422`；
* 返回单元按 **3′端阈值、连续互补阈值升序**排列；批次提交时的原阈值组合
  若落在请求范围内，对应单元标记 `is_batch_threshold: true`（范围外则无
  标记单元）；批次不存在返回 `404`。

```bash
curl -s -X POST http://localhost:8000/api/batches/f98cdd27-.../threshold-map \
  -H 'Content-Type: application/json' \
  -d '{
    "continuous_start": 3,
    "continuous_end": 6,
    "three_prime_start": 2,
    "three_prime_end": 4
  }'
```

请求字段：

| 字段 | 说明 |
| --- | --- |
| `continuous_start` / `continuous_end` | 连续互补阈值起止（含端点，各 1～200，起点 ≤ 终点） |
| `three_prime_start` / `three_prime_end` | 3′端连续互补阈值起止（含端点，各 1～200，起点 ≤ 终点） |

返回示例（2 引物 3 配对批次，原阈值 ct=4、tt=3）：

```json
{
  "batch_id": "f98cdd27-...",
  "continuous_start": 3, "continuous_end": 6,
  "three_prime_start": 2, "three_prime_end": 4,
  "cell_count": 12,
  "pair_count": 3,
  "batch_continuous_threshold": 4,
  "batch_three_prime_threshold": 3,
  "cells": [
    {
      "continuous_threshold": 3, "three_prime_threshold": 2,
      "high_count": 1, "medium_count": 0, "low_count": 2,
      "is_batch_threshold": false
    },
    {
      "continuous_threshold": 4, "three_prime_threshold": 3,
      "high_count": 1, "medium_count": 0, "low_count": 2,
      "is_batch_threshold": true
    }
  ]
}
```

* 每个单元的 `high_count + medium_count + low_count` 恒等于 `pair_count`；
* `cells` 完整覆盖请求范围内的全部阈值组合（每组合恰好一个单元），
  按 `three_prime_threshold`、`continuous_threshold` 升序；
* 原阈值单元的三级计数与批次摘要中的
  `high_risk_count` / `medium_risk_count` / `low_risk_count` 一致。

### `GET /api/batches/{base_id}/compare/{target_id}` — 批次风险对比

按引物**名称**构造稳定配对键，对比 base → target 两个已保存批次的配对风险：

* **配对键**：自二聚体保留单名（如 `F`）；异二聚体忽略批次内提交顺序，
  两名按字典序排列（如 `F|R`，与提交时谁先谁后无关）；
* **风险等级全序** `low < medium < high`：两批次共有的配对，等级上升进入
  `escalated`、下降进入 `deescalated`、不变不计入明细（计入
  `summary.unchanged_pair_count`）；仅存在于 target 的配对进入 `introduced`，
  仅存在于 base 的进入 `resolved`；
* **同名引物序列变化**仍按同一配对键比较，并在明细中标记
  `sequence_changed: true`（配对中任一同名引物序列不同即为 true）；
* 四类明细各自按配对名确定性排序，相同输入结果稳定可复现；
  任一批次不存在返回 `404`。

```bash
curl -s http://localhost:8000/api/batches/f98cdd27-.../compare/3fa9c1e2-...
```

返回示例：

```json
{
  "summary": {
    "base_batch_id": "f98cdd27-...",
    "target_batch_id": "3fa9c1e2-...",
    "base_pair_count": 15,
    "target_pair_count": 15,
    "common_pair_count": 10,
    "unchanged_pair_count": 7,
    "escalated_count": 1,
    "deescalated_count": 2,
    "introduced_count": 5,
    "resolved_count": 5
  },
  "escalated": [
    {
      "pair_key": "A1|A2",
      "name_a": "A1", "name_b": "A2", "is_self_dimer": false,
      "sequence_changed": true,
      "base_risk_level": "medium", "target_risk_level": "high"
    }
  ],
  "deescalated": [
    {
      "pair_key": "D1|D2",
      "name_a": "D1", "name_b": "D2", "is_self_dimer": false,
      "sequence_changed": true,
      "base_risk_level": "high", "target_risk_level": "low"
    }
  ],
  "introduced": [
    {
      "pair_key": "New",
      "name_a": "New", "name_b": "New", "is_self_dimer": true,
      "sequence_changed": false,
      "base_risk_level": null, "target_risk_level": "low"
    }
  ],
  "resolved": [
    {
      "pair_key": "A1|Gone",
      "name_a": "A1", "name_b": "Gone", "is_self_dimer": false,
      "sequence_changed": false,
      "base_risk_level": "high", "target_risk_level": null
    }
  ]
}
```

* `introduced` 明细的 `base_risk_level` 为 `null`，`resolved` 明细的
  `target_risk_level` 为 `null`；
* `summary.common_pair_count = unchanged + escalated + deescalated`，
  `introduced_count = target_pair_count - common_pair_count`，
  `resolved_count = base_pair_count - common_pair_count`。

## 配对结果字段含义

| 字段 | 含义 |
| --- | --- |
| `is_self_dimer` | 是否为自二聚体（同一引物与自身对齐） |
| `alignments_checked` | 实际枚举的错位数量，恒为 `len(top)+len(bottom)-1` |
| `best_shift` | 最危险对齐的带符号错位（见下） |
| `overlap_len` | 该对齐重叠碱基数 |
| `continuous_complementary_length` | 重叠区内最长连续互补长度 |
| `top3_run_length` / `bottom3_run_length` | 触及对应 3′端的连续互补长度（端未参与则为 0） |
| `three_prime_run_length` | 两端的最大值（`max(top3, bottom3)`） |
| `involved_three_prime_ends` | 3′端连续互补长度 > 0 的端名称 |
| `risk_level` | `high`（3′端达阈值）/ `medium`（仅连续互补达阈值）/ `low` |
| `risk_reasons` | 中文风险原因列表 |
| `alignment` | 等宽三行可视化；中间行 `\|`=互补，`.`=不互补，空格=未重叠 |

## 错位 `shift` 定义与同分规则

两条引物按反向平行对齐，top 固定 5′→3′，bottom 按 3′→5′ 打印。
`shift`（记为 d）是 **bottom 3′端碱基所在的全局列坐标**，top 的 5′端 T0 位于列 0：

下列三行与 API 的 `alignment` 输出格式一致（等宽，空格为未重叠，`|` 为互补）：

```
d = 0（齐平）        d > 0（bottom 3′伸入 top）   d < 0（bottom 3′悬空于左侧）
T0 T1 T2             T0 T1 T2                        T0 T1 T2
|  |  |              .  |  |  |                  |  |  .
   B2 B1 B0             B2 B1 B0              B2 B1 B0
（渲染为连续字符行，例如 d=0 时 top="TTT"、bottom="BBB"）
```

* `d = 0`：bottom 3′端碱基与 top 5′端 T0 同列；
* `d > 0`：bottom 3′端向 top 3′方向伸入 top 内部；
* `d < 0`：bottom 3′端悬空于 top 5′端左侧。

枚举范围为 `-(b-1) <= d <= a-1`，共 `a+b-1` 个对齐，全部检查，不做启发式剪枝。

**最危险对齐**按以下键取最大者：

1. `max(top3_run, bottom3_run)`（聚合酶延伸风险，最重要）；
2. 最长连续互补长度（双链稳定性）；
3. 重叠碱基数。

**同分固定偏移规则**（保证结果唯一、可复现）：风险键完全相同时依次取

1. `|d|` 最小者（悬空端最短、最稳定构象）；
2. `d` 为正者优先，其次 0，最后负值；
3. 再按 `d` 升序兜底（全序，无并列残留）。

## IUPAC 简并碱基互补规则

每个简并码映射为碱基集合（如 `R={A,G}`、`Y={C,T}`、`N={A,C,G,T}`）。
同列两碱基“可能互补”当且仅当一个碱基集合的 Watson-Crick 互补集合与另一个
集合**交集非空**。因此 `R` 与 `Y` 互补、`N` 与任意碱基互补、`S={G,C}` 与
`W={A,T}` 不互补。只要存在至少一种简并实例化方式可形成 WC 配对即判为互补。

---

## 本地开发（不使用 Docker）

```bash
pip install -r requirements.txt
export PRIMER_DB_PATH=./data/primer_dimer.db
python -m app.cli initdb
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install pytest
python -m pytest -q
```

## 项目结构

```
app/
  config.py     # 环境变量配置
  iupac.py      # IUPAC 碱基集合与互补矩阵
  dimer.py      # 错位枚举、连续互补/3′端长度、风险判定、同分规则
  pooling.py    # 异二聚体冲突图构建与 DSatur 风格确定性贪心分池
  augmentation.py # 候选增量冲突图上的确定性分支限界（容量受限最大独立集）
  ranking.py    # 替换候选的非支配分层（Pareto front，越小越优）
  threshold_map.py  # 阈值扫描地图的二维差分（imos）区间累加聚合
  compare.py    # 批次间稳定配对键构造与风险升降分类（escalated 等四类）
  schemas.py    # Pydantic 模型与输入校验
  db.py         # SQLite 建表、写入、查询
  services.py   # 批次计算编排（N 自二聚体 + N(N-1)/2 异二聚体）、分池、替换排名、候选增量筛选、阈值地图与对比编排
  main.py       # FastAPI 路由
  cli.py        # initdb 命令
tests/          # pytest 算法与 API 测试
Dockerfile
docker-compose.yml
docker-entrypoint.sh
requirements.txt
```
