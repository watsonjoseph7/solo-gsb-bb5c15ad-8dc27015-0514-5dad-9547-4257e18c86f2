"""引物二聚体对齐枚举与风险判定核心算法。

几何模型
--------
两条引物均按 5'→3' 给出；二聚体为反向平行（antiparallel）对齐，
可视化时 bottom 按 3'→5' 打印：

    top    5'  T0 T1  ... Ta-1 3'
                |  |
    bottom 3'  Bb-1 ... B1 B0  5'

带符号错位 ``shift``（记为 d）定义为 bottom 3′端碱基 B(b-1) 所在的全局列坐标，
top 的 T0 固定位于列 0：

* ``d = 0``：bottom 3′端与 top 5′端齐平；
* ``d > 0``：bottom 3′端伸入 top（向 top 3′方向）；
* ``d < 0``：bottom 3′端悬空于 top 5′端左侧。

配对关系（同列配对）：``t + b_idx = d + b - 1``。
全部存在重叠的错位范围为 ``-(b-1) <= d <= a-1``，共 ``a+b-1`` 个，
逐一枚举，不做采样。
"""

from __future__ import annotations

from dataclasses import dataclass

from .iupac import are_compatible


@dataclass(frozen=True)
class Alignment:
    """单个错位对齐的枚举结果。"""

    shift: int                      # 带符号错位 d（bottom 3′端的全局列坐标）
    overlap_len: int                # 重叠碱基数
    max_complementary_run: int      # 重叠区内最长连续互补长度
    top3_run: int                   # 触及 top 3′端的连续互补长度（无则 0）
    bottom3_run: int                # 触及 bottom 3′端的连续互补长度（无则 0）
    top_line: str                   # top 可视化行（5'→3'，等宽）
    match_line: str                 # 中间行：| 可能互补，. 不互补
    bottom_line: str                # bottom 反向平行展示行（左→右为 3'→5'）


@dataclass(frozen=True)
class PairResult:
    """一对引物全部错位枚举后的最危险对齐与风险结论。"""

    name_top: str
    name_bottom: str
    is_self_dimer: bool
    alignments_checked: int
    best_shift: int
    overlap_len: int
    max_complementary_run: int
    top3_run: int
    bottom3_run: int
    terminal_run: int               # max(top3_run, bottom3_run)
    involved_three_prime_ends: tuple[str, ...]
    risk_level: str                 # high / medium / low
    risk_reasons: tuple[str, ...]
    top_line: str
    match_line: str
    bottom_line: str


def _enumerate_alignment(top: str, bottom: str, d: int) -> Alignment:
    """枚举给定错位 d 下的一条反向平行对齐。"""
    a, b = len(top), len(bottom)

    # 同列配对 t + b_idx = d + b - 1：
    #   b_idx <= b-1  -> t >= d
    #   b_idx >= 0    -> t <= d + b - 1
    t_lo = max(0, d)
    t_hi = min(a - 1, d + b - 1)

    pairs_top: list[str] = []
    pairs_bottom: list[str] = []
    matches: list[bool] = []
    for t in range(t_lo, t_hi + 1):
        b_idx = d + b - 1 - t
        pairs_top.append(top[t])
        pairs_bottom.append(bottom[b_idx])
        matches.append(are_compatible(top[t], bottom[b_idx]))

    # 最长连续互补
    best_run = run = 0
    for ok in matches:
        run = run + 1 if ok else 0
        if run > best_run:
            best_run = run

    # bottom 3′端（B(b-1)）在 d>=0 时出现于重叠区第一列（t=d 为最小重叠 t）
    bottom3_run = 0
    if d >= 0 and matches and matches[0]:
        for ok in matches:
            if ok:
                bottom3_run += 1
            else:
                break

    # top 3′端（T(a-1)）在 d >= a-b 时出现于重叠区最后一列
    top3_run = 0
    if d >= a - b and matches and matches[-1]:
        for ok in reversed(matches):
            if ok:
                top3_run += 1
            else:
                break

    # 等宽可视化。全局列区间 [x_lo, x_hi)。
    # top 占 [0, a)；bottom 占 [d, d+b)。
    x_lo = min(0, d)
    x_hi = max(a, d + b)
    width = x_hi - x_lo

    top_cells = [" "] * width
    bottom_cells = [" "] * width
    match_cells = [" "] * width

    for t in range(a):
        top_cells[t - x_lo] = top[t]
    # B(b_idx) 位于列 d+b-1-b_idx；打印左→右为 B(b-1)..B(0)，即 3'→5'
    for b_idx in range(b):
        bottom_cells[d + b - 1 - b_idx - x_lo] = bottom[b_idx]
    for idx, ok in enumerate(matches):
        col = (t_lo + idx) - x_lo
        match_cells[col] = "|" if ok else "."

    return Alignment(
        shift=d,
        overlap_len=len(matches),
        max_complementary_run=best_run,
        top3_run=top3_run,
        bottom3_run=bottom3_run,
        top_line="".join(top_cells),
        match_line="".join(match_cells),
        bottom_line="".join(bottom_cells),
    )


def _tie_key(aln: Alignment) -> tuple:
    """同分固定偏移规则（值越小越优先）。

    风险键完全相同的同分对齐，按以下顺序取唯一结果：
    1. ``|shift|`` 最小（错位最少、双链悬空端最短）；
    2. shift 为正者优先（+d 表示 bottom 3′端伸入 top，是更常规的延伸构象），
       其次 0，最后负值；
    3. 仍相同（自二聚体的镜像情形不会走到这里）则按 shift 升序兜底，保证全序。
    """
    sign = (aln.shift > 0) - (aln.shift < 0)  # 正:1 零:0 负:-1
    return (abs(aln.shift), -sign, aln.shift)


def _risk_key(aln: Alignment) -> tuple:
    """“最危险对齐”排序键（越大越危险）。

    1. 任一 3′端连续互补长度（聚合酶可沿二聚体延伸，首要风险）；
    2. 最长连续互补长度（双链热稳定性）；
    3. 重叠碱基总数（互补机会总量）。
    """
    return (
        max(aln.top3_run, aln.bottom3_run),
        aln.max_complementary_run,
        aln.overlap_len,
    )


def _classify_risk(
    aln: Alignment, continuous_threshold: int, three_prime_threshold: int
) -> tuple[str, list[str]]:
    terminal = max(aln.top3_run, aln.bottom3_run)
    reasons: list[str] = []

    if terminal >= three_prime_threshold:
        ends: list[str] = []
        if aln.top3_run >= three_prime_threshold:
            ends.append(f"top 3′端连续互补 {aln.top3_run} bp")
        if aln.bottom3_run >= three_prime_threshold:
            ends.append(f"bottom 3′端连续互补 {aln.bottom3_run} bp")
        reasons.append(
            "3′端连续互补达到阈值，存在聚合酶延伸风险：" + "，".join(ends)
        )
    if aln.max_complementary_run >= continuous_threshold:
        reasons.append(
            f"最长连续互补 {aln.max_complementary_run} bp 达到阈值 "
            f"{continuous_threshold} bp，双链结合稳定"
        )

    if terminal >= three_prime_threshold:
        level = "high"
    elif aln.max_complementary_run >= continuous_threshold:
        level = "medium"
    else:
        level = "low"
        reasons.append(
            f"最长连续互补 {aln.max_complementary_run} bp、3′端连续 {terminal} bp "
            "均未达到阈值，风险较低"
        )
    return level, reasons


def analyze_pair(
    top_name: str,
    top: str,
    bottom_name: str,
    bottom: str,
    continuous_threshold: int,
    three_prime_threshold: int,
    is_self_dimer: bool = False,
) -> PairResult:
    """枚举一对引物全部 ``a+b-1`` 个错位对齐，返回最危险对齐。"""
    a, b = len(top), len(bottom)
    alignments = [
        _enumerate_alignment(top, bottom, d)
        for d in range(-(b - 1), a)
    ]

    best = max(
        alignments,
        # 风险键越大越危险；同分时固定偏移规则（tie 键越小越优先 -> 取负）
        key=lambda aln: (_risk_key(aln), tuple(-x for x in _tie_key(aln))),
    )

    level, reasons = _classify_risk(
        best, continuous_threshold, three_prime_threshold
    )

    ends: list[str] = []
    if best.top3_run > 0:
        ends.append(f"{top_name}:3′")
    if best.bottom3_run > 0:
        ends.append(f"{bottom_name}:3′")

    return PairResult(
        name_top=top_name,
        name_bottom=bottom_name,
        is_self_dimer=is_self_dimer,
        alignments_checked=len(alignments),
        best_shift=best.shift,
        overlap_len=best.overlap_len,
        max_complementary_run=best.max_complementary_run,
        top3_run=best.top3_run,
        bottom3_run=best.bottom3_run,
        terminal_run=max(best.top3_run, best.bottom3_run),
        involved_three_prime_ends=tuple(ends),
        risk_level=level,
        risk_reasons=tuple(reasons),
        top_line=best.top_line,
        match_line=best.match_line,
        bottom_line=best.bottom_line,
    )
