"""批次内冲突修复的最小代价移除集求解（确定性分支限界）。

冲突图
------
顶点为请求中按批次位置升序给出的 2～28 条批次内引物（顶点号即该引物在
本模块输入 ``edges`` 中的位置下标）。命中所选风险等级的配对产生一条冲突边：

* 引物的**自二聚体**命中所选风险等级时为**自环** (v, v)，自环只能靠移除
  该引物自身覆盖，因此自环节点被强制移除；若该节点同时受保护则冲突不可消除；
* 两条引物之间的命中配对记普通无向边 (u, v)；
* 两个**受保护**节点之间的普通边无法被任何移除集覆盖（两端都不可移除），
  同样属于不可消除冲突。

问题模型
--------
在“自环节点必移除、受保护节点不可移除”的约束下，求覆盖全部冲突边的
**最小权顶点覆盖**：每个顶点带正整数移除代价。候选移除集依次按
① **移除总代价**（小者优先）、② **移除数量**（少者优先）、
③ **批次位置序列**（即移除下标元组，字典序小者优先）比较，故全序无并列。

分支限界
--------
自环与“受保护端的唯一可选端”先做确定性强制传播（含级联），只在剩余
未决顶点上按下标升序做“**移除优先 / 保留其后**”的深度优先枚举。
“移除优先 + 下标升序”使解按移除位置元组的字典序被访问，故在同代价、
同数量的最优解中字典序最小者为唯一返回。剪枝手段：

* 下界：剩余边的贪心匹配上，每条匹配边至少贡献其两端的最小代价；
  移除数量下界为匹配边数。代价或数量不可能优于当前最优即剪枝；
* 保留一个顶点时，其全部剩余邻点被强制移除，遇已保留邻点立即回溯。

输入规模 2～28，最坏情形亦可瞬时完成。
"""

from __future__ import annotations

from typing import NamedTuple


class InfeasibleRepair(Exception):
    """保护约束使问题无解：存在任何移除集都无法覆盖的冲突边。"""

    def __init__(self, edges: list[tuple[int, int]]):
        # edges 为不可消除冲突边 (u, v)，u <= v；u == v 即自环
        self.edges = sorted(set(edges))
        super().__init__(f"不可消除冲突边 {len(self.edges)} 条")


class RepairPlan(NamedTuple):
    """精确求解结果：移除顶点下标（升序）与保留顶点下标（升序）。"""

    removed: tuple[int, ...]
    retained: tuple[int, ...]


def minimum_cost_removal(
    edges: list[tuple[int, int]],
    *,
    costs: list[int],
    protected: frozenset[int],
) -> RepairPlan:
    """求覆盖 ``edges`` 全部冲突边的最小代价移除集。

    ``edges`` 允许自环 (v, v)；普通边端点顺序任意，模块内部按 (u, v)
    u < v 规范化。顶点编号为 0..n-1，``costs[v]`` 为正整数移除代价，
    ``protected`` 为不可移除顶点。保护约束令问题无解时抛出
    :class:`InfeasibleRepair`，并在 ``edges`` 中列出全部不可消除冲突
    （受保护自环、两端均受保护的普通边）。
    """
    n = len(costs)

    # 自环去重；普通边规范化为 u < v
    loops: set[int] = set()
    undirected: set[tuple[int, int]] = set()
    for raw_u, raw_v in edges:
        u, v = (raw_u, raw_v) if raw_u <= raw_v else (raw_v, raw_u)
        if u == v:
            loops.add(u)
        else:
            undirected.add((u, v))

    # 不可消除冲突：受保护自环、两端均受保护的普通边
    impossible = sorted(
        {(v, v) for v in loops if v in protected}
        | {(u, v) for u, v in undirected if u in protected and v in protected}
    )
    if impossible:
        raise InfeasibleRepair(impossible)

    # 确定性强制传播：自环节点必移除；保留节点（受保护节点，以及后续级联中
    # 被保留的节点）的邻点必移除。级联只可能因“两个保留节点之间有边”而失败，
    # 但初始保留集合 protected 之间的边已在上面判为不可消除
    removed: set[int] = {v for v in loops if v not in protected}
    kept: set[int] = set(protected)

    changed = True
    while changed:
        changed = False
        for u, v in undirected:
            u_removed, v_removed = u in removed, v in removed
            if u_removed or v_removed:
                continue  # 边已被覆盖
            u_kept, v_kept = u in kept, v in kept
            if u_kept and v_kept:
                raise InfeasibleRepair([(u, v)])
            if u_kept:
                removed.add(v)
                changed = True
            elif v_kept:
                removed.add(u)
                changed = True

    # 剩余未决顶点（自环节点与受保护节点均已离开未决集合）
    residual_adjacency: list[set[int]] = [set() for _ in range(n)]
    for u, v in undirected:
        if u in removed or v in removed or u in kept or v in kept:
            continue
        residual_adjacency[u].add(v)
        residual_adjacency[v].add(u)
    free = sorted(v for v in range(n) if residual_adjacency[v])

    forced_cost = sum(costs[v] for v in removed)

    # 贪心初始上界：反复取剩余代价最低的顶点覆盖一条未覆盖边，
    # 得到一个可行移除集（加上强制移除），及早建立剪枝门槛
    greedy: set[int] = set()
    pending = {(u, v) for u, v in undirected
               if u not in removed and v not in removed}
    while pending:
        u, v = min(pending)
        pick = u if (costs[u], u) <= (costs[v], v) else v
        greedy.add(pick)
        pending = {e for e in pending if pick not in e}

    best_removed: set[int] = removed | greedy
    best_key: tuple[int, int, tuple[int, ...]] = (
        forced_cost + sum(costs[v] for v in greedy),
        len(best_removed),
        tuple(sorted(best_removed)),
    )

    chosen: list[int] = []  # 分支中移除的 free 顶点（按下标升序追加）

    def lower_bound(current_cost: int) -> tuple[int, int]:
        """贪心匹配给出的 (总代价, 移除数量) 下界。

        每条匹配边的两个端点在最终覆盖中至少移除一个，且该端点代价不低于
        两端代价的较小值；匹配边两两不共端点，故各边的最小端点代价可相加，
        匹配边数同时是额外移除数量的下界。
        """
        matched: set[int] = set()
        bound_cost = current_cost
        bound_count = len(removed) + len(chosen)
        for v in free:
            if v in matched or v in kept or v in chosen:
                continue
            for u in sorted(residual_adjacency[v]):
                if u in matched or u in kept or u in chosen:
                    continue
                # (v, u) 两个端点均未决且未在本次分支移除
                matched.add(v)
                matched.add(u)
                bound_cost += min(costs[v], costs[u])
                bound_count += 1
                break
        return bound_cost, bound_count

    def search(i: int) -> None:
        nonlocal best_key, best_removed

        if i == len(free):
            removed_tuple = tuple(sorted(removed | set(chosen)))
            key = (
                forced_cost + sum(costs[v] for v in chosen),
                len(removed_tuple),
                removed_tuple,
            )
            if key < best_key:
                best_key = key
                best_removed = set(removed_tuple)
            return

        v = free[i]
        if v in chosen or v in kept:
            search(i + 1)
            return

        # 下界剪枝（代价、数量任一维不可能优于当前最优即放弃整棵子树）
        current_cost = forced_cost + sum(costs[x] for x in chosen)
        bound_cost, bound_count = lower_bound(current_cost)
        if bound_cost > best_key[0] or (
            bound_cost == best_key[0] and bound_count > best_key[1]
        ):
            return

        neighbors = residual_adjacency[v]

        # 分支一：移除 v（移除优先：同代价同数量时先访问字典序更小的位置序列）
        chosen.append(v)
        search(i + 1)
        chosen.pop()

        # 分支二：保留 v，其全部剩余邻点必须移除（强制传播）；
        # 邻点中已有保留节点则该边无法覆盖，立即回溯
        if any(u in kept for u in neighbors):
            return
        forced_here = [u for u in sorted(neighbors) if u not in chosen]
        chosen.extend(forced_here)
        kept.add(v)
        search(i + 1)
        kept.discard(v)
        for u in forced_here:
            chosen.remove(u)

    search(0)

    removed_final = tuple(sorted(best_removed))
    retained_final = tuple(v for v in range(n) if v not in best_removed)
    return RepairPlan(removed=removed_final, retained=retained_final)
