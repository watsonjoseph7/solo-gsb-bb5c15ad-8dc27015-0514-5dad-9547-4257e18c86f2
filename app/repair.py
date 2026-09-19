"""批次内引物冲突修复：含自环冲突图上的最小代价移除集（确定性分支限界）。

冲突图
------
顶点为请求指定的批次内引物（顶点号即按批次位置升序重排后的序号）。所选
风险等级的配对产生冲突：

* 自二聚体达到所选风险等级 → 该顶点带**自环**，自环节点必须被移除
  （任何移除集都不可能“覆盖”自环以外的出路）；
* 两引物间异二聚体达到所选风险等级 → 普通无向边，移除任一端即可覆盖。

受保护顶点不可移除：自环落在受保护顶点、或冲突边两端均受保护时问题无解
（服务层据此返回 409）；恰好一端受保护的边等价于另一端被强制移除。

分支限界
--------
预处理把自环节点与受保护顶点的邻居计入**强制移除**后，在剩余无向图上求
最小代价顶点覆盖，依次以 ① 移除总成本、② 移除数量、③ 顶点号（批次位置）
升序元组 决胜——强制部分在所有可行解中相同，故只需对分支限界部分比较该
三元组。搜索按顶点号升序做“包含优先 / 排除其后”的深度优先枚举：排除顶点
时其在剩余图中的全部邻居被强制包含。每个搜索节点用贪心团覆盖（补图着色）
计算剩余图的加权顶点覆盖下界——每个团至少贡献“团内代价和 − 团内最大代价”
的成本与“团大小 − 1”的顶点数——与当前最优的 (成本, 数量) 前缀比较剪枝；
当成本与数量下界均与当前最优打平时，再用“chosen 并入剩余活跃顶点中最小
的 k 个”的乐观元组估计做决胜剪枝。叶子处以完整三元组决胜，保证全局最优
且决胜规则无并列残留。输入规模 2～28，最坏情形亦可瞬时完成。
"""

from __future__ import annotations

from typing import NamedTuple


class VertexCover(NamedTuple):
    """精确求解结果：全部移除顶点（含强制，升序）与其中被强制移除的顶点。"""

    removed: tuple[int, ...]
    forced: frozenset[int]


def _clique_cover_bounds(
    adjacency: list[set[int]], costs: list[int], vertices: frozenset[int]
) -> tuple[int, int]:
    """贪心团覆盖给出的加权顶点覆盖 (成本下界, 顶点数下界)。

    按顶点号升序贪心地把每个顶点放入第一个“与团内全体顶点相邻”的类，
    放不下则新开一类（等价于补图贪心着色）。大小为 k 的团在任何顶点覆盖
    中至少选中 k-1 个顶点，故该类至少贡献“团内代价和 − 团内最大代价”的
    成本与 k-1 个顶点；团之间顶点不相交，各项下界可直接相加。
    """
    classes: list[list[int]] = []
    for v in sorted(vertices):
        for cls in classes:
            if all(u in adjacency[v] for u in cls):
                cls.append(v)
                break
        else:
            classes.append([v])
    cost_lb = 0
    count_lb = 0
    for cls in classes:
        if len(cls) > 1:
            cost_lb += sum(costs[u] for u in cls) - max(costs[u] for u in cls)
            count_lb += len(cls) - 1
    return cost_lb, count_lb


def minimum_cost_cover(
    *,
    costs: list[int],
    protected: frozenset[int],
    self_loops: frozenset[int],
    edges: list[tuple[int, int]],
) -> VertexCover | None:
    """确定性分支限界求最小代价顶点覆盖（含自环与受保护约束）。

    ``costs[v]`` 为顶点 v 的正整数移除代价；``protected`` 为不可移除顶点；
    ``self_loops`` 为带自环顶点（必须移除）；``edges`` 为普通无向边
    ``(u, v)``（u < v）。受保护约束使问题无解时返回 None（自环落在受保护
    顶点，或某条边两端均受保护）；否则返回唯一最优解：依次以移除总成本、
    移除数量、顶点号升序元组决胜。
    """
    n = len(costs)
    adjacency: list[set[int]] = [set() for _ in range(n)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)

    # 可行性：受保护顶点无法参与任何移除集
    if self_loops & protected:
        return None
    if any(u in protected and v in protected for u, v in edges):
        return None

    # 强制移除：自环节点，以及受保护顶点在冲突边另一端的邻居
    forced: set[int] = set(self_loops)
    for u, v in edges:
        if u in protected:
            forced.add(v)
        if v in protected:
            forced.add(u)
    forced -= protected

    initial = frozenset(set(range(n)) - protected - forced)

    best_key: tuple[int, int, tuple[int, ...]] | None = None
    best: tuple[int, ...] = ()
    chosen: list[int] = []
    chosen_cost = 0

    def search(remaining: frozenset[int]) -> None:
        nonlocal best_key, best, chosen_cost

        # 活跃顶点：剩余图中仍有关联冲突边者（升序）；无关联边的顶点默认
        # 排除（正代价下不可能改善任何目标），也不会再被选中
        active = [v for v in sorted(remaining) if adjacency[v] & remaining]

        if not active:
            key = (chosen_cost, len(chosen), tuple(sorted(chosen)))
            if best_key is None or key < best_key:
                best_key = key
                best = tuple(sorted(chosen))
            return

        # 下界剪枝：仅当 (成本, 数量) 下界严格劣于当前最优时才可剪，
        # 相等时子树仍可能藏有顶点号元组更小的同优解
        cost_lb, count_lb = _clique_cover_bounds(adjacency, costs, remaining)
        if best_key is not None:
            cost_lo = chosen_cost + cost_lb
            count_lo = len(chosen) + count_lb
            if cost_lo > best_key[0]:
                return
            if cost_lo == best_key[0]:
                if count_lo > best_key[1]:
                    return
                if count_lo == best_key[1]:
                    # 成本与数量都只能恰好打平，决胜落在顶点号元组上。
                    # 本子树可达的最小元组 = chosen 并入活跃顶点中最小的
                    # k 个（k 为达到最优数量所需的顶点数，k <= 活跃顶点数）；
                    # 若该乐观估计已不优于当前最优，整棵子树不可能更优
                    k = best_key[1] - len(chosen)
                    optimistic = tuple(sorted(chosen + active[:k]))
                    if optimistic >= best_key[2]:
                        return

        v = active[0]
        rest = remaining - {v}

        # 分支一：包含 v，其关联边全部被覆盖
        chosen.append(v)
        chosen_cost += costs[v]
        search(rest)
        chosen_cost -= costs[v]
        chosen.pop()

        # 分支二：排除 v，其在剩余图中的全部邻居被强制包含
        neighbors = adjacency[v] & rest
        added = sorted(neighbors)
        chosen.extend(added)
        chosen_cost += sum(costs[u] for u in added)
        search(rest - neighbors)
        chosen_cost -= sum(costs[u] for u in added)
        del chosen[len(chosen) - len(added):]

    search(initial)
    removed = tuple(sorted(set(best) | forced))
    return VertexCover(removed=removed, forced=frozenset(forced))
