"""异二聚体冲突图上的引物分池（DSatur 风格确定性贪心）。

冲突图
------
顶点为批次内引物（顶点号即批次位置），边为需要隔离风险等级的**异二聚体**
配对（自二聚体在进入本模块前已被忽略）。两端点落在同一池中的边即
“同池冲突”。

贪心规则（全部比较均为全序，相同输入必得相同结果）
--------------------------------------------------
1. 选择下一个引物（未分配顶点中取排序键最大者）：
   a. **已见邻池数**（saturation degree）：相邻引物已占据的不同池的数量；
   b. **冲突度**（degree）：该引物在冲突图中的边数；
   c. **批次位置**：位置靠前者优先（稳定兜底）。
2. 确定放置池（仅考虑未满池，取排序键最小者）：
   a. **新增冲突数**：放入该池后新产生的同池冲突边数；
   b. **池内人数**：池内已有引物数，少者优先；
   c. **池号**：小者优先。

调用方须保证 ``pool_count * pool_capacity >= 引物数``；该条件下贪心每一步
都至少存在一个未满池，全部引物必被分配。
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ConflictEdge:
    """一条同池冲突边（对应一对异二聚体）。"""

    position_a: int
    position_b: int
    risk_level: str


@dataclass(frozen=True)
class PoolPlan:
    """分池结果。

    assignments: 顶点（批次位置）→ 池号（0 起）；
    pool_members: 池号 → 成员位置列表（按批次位置升序）；
    conflicts: (池号, 冲突边) 列表，按池号、两端位置升序。
    """

    assignments: dict[int, int]
    pool_members: list[list[int]]
    conflicts: list[tuple[int, ConflictEdge]]

    @property
    def total_conflicts(self) -> int:
        return len(self.conflicts)


def build_pool_plan(
    *,
    primer_count: int,
    edges: list[tuple[int, int, str]],
    pool_count: int,
    pool_capacity: int,
) -> PoolPlan:
    """在异二聚体冲突图上执行 DSatur 风格贪心分池。

    edges 元素为 ``(位置u, 位置v, 风险等级)``，u/ v 为批次位置，无重复无向边。
    """
    n = primer_count
    adjacency: list[set[int]] = [set() for _ in range(n)]
    edge_risk: dict[tuple[int, int], str] = {}
    for u, v, level in edges:
        if u == v:
            # 自二聚体不进入冲突图（双保险，正常应由上层过滤）
            continue
        a, b = (u, v) if u < v else (v, u)
        adjacency[a].add(b)
        adjacency[b].add(a)
        edge_risk[(a, b)] = level

    # color[v] = 池号；-1 表示未分配
    color = [-1] * n
    sizes = [0] * pool_count
    # 已见邻池数与邻居已占池集合（增量维护，避免每轮重扫）
    saturation = [0] * n
    neighbor_pools: list[set[int]] = [set() for _ in range(n)]
    uncolored = set(range(n))

    while uncolored:
        # 选择：已见邻池数 ↓、冲突度 ↓、批次位置 ↑
        v = max(
            uncolored,
            key=lambda x: (saturation[x], len(adjacency[x]), -x),
        )

        # 该引物放入各池的新增冲突数 = 已在该池中的相邻引物数
        added_conflicts = [0] * pool_count
        for u in adjacency[v]:
            c = color[u]
            if c >= 0:
                added_conflicts[c] += 1

        # 放置：新增冲突数 ↑少、池内人数 ↑少、池号 ↑小（均取最小）
        chosen = min(
            (c for c in range(pool_count) if sizes[c] < pool_capacity),
            key=lambda c: (added_conflicts[c], sizes[c], c),
        )

        color[v] = chosen
        sizes[chosen] += 1
        uncolored.remove(v)
        for u in adjacency[v]:
            if color[u] < 0 and chosen not in neighbor_pools[u]:
                neighbor_pools[u].add(chosen)
                saturation[u] += 1

    pool_members: list[list[int]] = [[] for _ in range(pool_count)]
    for position, c in enumerate(color):
        pool_members[c].append(position)

    conflicts: list[tuple[int, ConflictEdge]] = []
    for (a, b), level in edge_risk.items():
        if color[a] == color[b]:
            conflicts.append((color[a], ConflictEdge(a, b, level)))
    conflicts.sort(key=lambda item: (item[0], item[1].position_a, item[1].position_b))

    return PoolPlan(
        assignments={position: c for position, c in enumerate(color)},
        pool_members=pool_members,
        conflicts=conflicts,
    )
