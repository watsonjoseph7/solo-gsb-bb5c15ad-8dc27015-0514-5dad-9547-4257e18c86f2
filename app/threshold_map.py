"""阈值扫描地图：二维差分（imos）区间累加聚合。

对已保存批次的每个配对，只取两个已落库的标量指标——最长连续互补长度
``c`` 与 3′端连续互补长度 ``t``——在请求给定的阈值网格上按现有判定规则
统计每个阈值组合下的三级配对数量：

* ``t >= 3′端阈值``         → high（聚合酶延伸风险，最优先）；
* 否则 ``c >= 连续互补阈值`` → medium（双链结合稳定）；
* 否则                      → low。

朴素做法需对每个配对遍历全部网格单元（O(配对数 × 单元数)）。本模块改用
**二维差分数组**：每个配对对 high、medium 计数各贡献一个轴对齐矩形阈值
区域，只在矩形四个角点 O(1) 打差分标记，最后用一次二维前缀和物化整张
地图，总复杂度 O(配对数 + 单元数)。
"""

from __future__ import annotations


def _add_rect(diff: list[list[int]], r0: int, c0: int, r1: int, c1: int) -> None:
    """闭区间矩形 [r0..r1] × [c0..c1] 整体 +1（仅打四个角点差分标记）。"""
    diff[r0][c0] += 1
    diff[r0][c1 + 1] -= 1
    diff[r1 + 1][c0] -= 1
    diff[r1 + 1][c1 + 1] += 1


def _materialize(diff: list[list[int]], rows: int, cols: int) -> list[list[int]]:
    """二维前缀和物化差分数组，返回 rows × cols 的计数网格。"""
    grid = [[0] * cols for _ in range(rows)]
    for r in range(rows):
        for c in range(cols):
            v = diff[r][c]
            if r:
                v += grid[r - 1][c]
            if c:
                v += grid[r][c - 1]
            if r and c:
                v -= grid[r - 1][c - 1]
            grid[r][c] = v
    return grid


def build_threshold_map(
    metrics: list[tuple[int, int]],
    *,
    continuous_start: int,
    continuous_end: int,
    three_prime_start: int,
    three_prime_end: int,
) -> list[list[tuple[int, int, int]]]:
    """在阈值网格上聚合三级配对数量。

    metrics 为每个配对的 ``(最长连续互补 c, 3′端连续互补 t)``。
    返回 ``grid[row][col] = (high, medium, low)``：行对应 3′端阈值
    ``three_prime_start + row``，列对应连续互补阈值
    ``continuous_start + col``，均按升序。每个单元的三级计数之和恒等于
    配对总数。
    """
    rows = three_prime_end - three_prime_start + 1
    cols = continuous_end - continuous_start + 1

    high_diff = [[0] * (cols + 1) for _ in range(rows + 1)]
    medium_diff = [[0] * (cols + 1) for _ in range(rows + 1)]

    for c, t in metrics:
        # high：3′端阈值 <= t 的全部单元（行前缀区间，列方向全覆盖）
        if t >= three_prime_start:
            _add_rect(
                high_diff, 0, 0,
                min(three_prime_end, t) - three_prime_start, cols - 1,
            )
        # medium：3′端阈值 > t 且连续互补阈值 <= c 的轴对齐矩形
        if c >= continuous_start and t < three_prime_end:
            _add_rect(
                medium_diff,
                max(three_prime_start, t + 1) - three_prime_start, 0,
                rows - 1, min(continuous_end, c) - continuous_start,
            )

    high_grid = _materialize(high_diff, rows, cols)
    medium_grid = _materialize(medium_diff, rows, cols)

    total = len(metrics)
    return [
        [
            (
                high_grid[r][c],
                medium_grid[r][c],
                total - high_grid[r][c] - medium_grid[r][c],
            )
            for c in range(cols)
        ]
        for r in range(rows)
    ]
