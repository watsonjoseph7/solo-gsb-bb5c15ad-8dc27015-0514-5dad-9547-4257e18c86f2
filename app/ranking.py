"""候选替换引物的非支配分层（Pareto front）。

每个候选带一个三维指标向量 ``(high 配对数, 最大 3′连续互补, 最大连续互补)``，
三个维度均为**越小越优**。候选 i 支配（dominate）候选 j 当且仅当 i 在所有
维度上不劣于 j，且至少一个维度严格更优。

分层采用快速非支配排序：第 1 层为全体非支配候选；移除第 1 层后，剩余候选
中的非支配集合为第 2 层，依此类推。同层候选之间互不支配，由调用方按输入
顺序稳定排序。输入规模 2～50，直接 O(n²) 两两比较即可。
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TypeVar

T = TypeVar("T")


def dominates(a: tuple[int, ...], b: tuple[int, ...]) -> bool:
    """a 支配 b：所有维度 a<=b 且至少一个维度 a<b（越小越优）。"""
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))


def non_dominated_layers(
    items: Sequence[T], metrics: Callable[[T], tuple[int, ...]]
) -> list[list[T]]:
    """把 items 划分为若干非支配层，返回逐层列表（每层保持 items 中的相对顺序）。

    metrics 把单个元素映射为同维整数指标向量（各维越小越优）。
    指标向量完全相同的元素互不支配，必然落在同一层。
    """
    n = len(items)
    vectors = [metrics(item) for item in items]

    # dominated_by[i]：被 i 支配的元素下标集合
    dominated_by: list[set[int]] = [set() for _ in range(n)]
    # active_dominators[j]：仍在剩余集合中、且支配 j 的元素数
    active_dominators = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and dominates(vectors[i], vectors[j]):
                dominated_by[i].add(j)
                active_dominators[j] += 1

    layers: list[list[T]] = []
    remaining = set(range(n))
    while remaining:
        # 剩余集合中没有支配者的元素构成下一层；按输入下标稳定输出
        front = [i for i in range(n) if i in remaining and active_dominators[i] == 0]
        layers.append([items[i] for i in front])
        for i in front:
            remaining.discard(i)
            for j in dominated_by[i]:
                active_dominators[j] -= 1
    return layers
