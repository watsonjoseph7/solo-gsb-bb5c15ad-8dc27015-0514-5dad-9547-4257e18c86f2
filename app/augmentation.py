"""候选增量的冲突图与确定性分支限界精确求解。

冲突图
------
顶点为请求中按顺序给出的候选引物（顶点号即候选下标，0 起）。一条 high 或
medium 配对产生一条冲突边：

* 候选自二聚体为 high/medium 时，该候选与自身冲突（自环），自环节点不可能
  被任何无冲突集合接纳；
* 候选与原批次引物的 high/medium 配对同样记为该候选的自环（原批次引物固定
  在图外且无条件“在场”，效果等价于候选被固定禁止）；
* 候选两两之间的 high/medium 配对记为普通无向边。

分支限界
--------
在上述冲突图上求**不超过 ``max_selected``** 的最大独立集（自环节点排除在
独立集之外）。顶点按下标升序做“包含优先 / 排除其后”的深度优先枚举，搜索
顺序本身即按选中下标元组的字典序访问各独立集，故首次达到最优规模的解就是
同规模解中的字典序最小者，全序无并列。每棵子树用“剩余容量、剩余可入选项数、
贪心团覆盖（补图着色）独立数上界”三者最小值剪枝；输入规模 2～18，最坏情形
亦可瞬时完成。
"""

from __future__ import annotations

from typing import NamedTuple


class IndependentSet(NamedTuple):
    """精确求解结果：选中顶点下标（升序）与因自环被排除的顶点下标（升序）。"""

    selected: tuple[int, ...]
    self_blocked: frozenset[int]


def _clique_cover_bound(adjacency: list[set[int]], vertices: frozenset[int]) -> int:
    """贪心团覆盖得到的独立数上界：同类顶点两两相邻（构成团）。

    独立集与每个团至多交于一个顶点，故团覆盖的团数是独立数的合法上界
    （等价于在补图上做贪心着色）。按顶点下标升序贪心地把每个顶点放入
    第一个“与团内全体顶点相邻”的类，放不下则新开一类，返回类数。
    """
    classes: list[list[int]] = []
    for v in sorted(vertices):
        for cls in classes:
            if all(u in adjacency[v] for u in cls):
                cls.append(v)
                break
        else:
            classes.append([v])
    return len(classes)


def maximum_independent_set(
    adjacency: list[set[int]], max_selected: int
) -> IndependentSet:
    """确定性分支限界求容量受限的最大独立集。

    ``adjacency[v]`` 为顶点 v 的全部邻居；邻居含 v 自身即自环（被原批次或
    自身自二聚体阻断）。同规模时选中下标元组字典序最小的解为唯一返回。
    """
    n = len(adjacency)
    self_blocked = frozenset(v for v in range(n) if v in adjacency[v])

    best: tuple[int, ...] = ()
    chosen: list[int] = []

    def search(remaining: frozenset[int]) -> None:
        nonlocal best

        if not remaining:
            if len(chosen) > len(best):
                best = tuple(chosen)
            return

        # 上界剪枝：剩余容量、剩余可入选项数、贪心团覆盖独立数上界三者取最小
        room = max_selected - len(chosen)
        if room <= 0:
            if len(chosen) > len(best):
                best = tuple(chosen)
            return
        if len(chosen) + min(
            room, len(remaining), _clique_cover_bound(adjacency, remaining)
        ) <= len(best):
            return

        # 按顶点下标升序分支：先“包含”后“排除”，搜索顺序即选中元组字典序，
        # 因此首次达到最优规模的解必为同规模解中字典序最小者
        v = min(remaining)
        rest = remaining - {v}

        # 分支一：包含 v（自环节点不会出现在 remaining 中），其邻居一并剔除
        chosen.append(v)
        search(rest - adjacency[v])
        chosen.pop()
        if len(best) == max_selected:
            return  # 已达容量上限，任何解都不可能更大

        # 分支二：排除 v
        search(rest)

    initial = frozenset(set(range(n)) - self_blocked)
    search(initial)
    return IndependentSet(selected=best, self_blocked=self_blocked)
