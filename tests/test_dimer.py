"""核心算法与 API 的测试。运行：python -m pytest -q"""

from __future__ import annotations

import pytest

from app.dimer import _enumerate_alignment, analyze_pair
from app.iupac import are_compatible, validate_iupac


# ---------- IUPAC ----------

def test_standard_base_complement():
    assert are_compatible("A", "T")
    assert are_compatible("T", "A")
    assert are_compatible("G", "C")
    assert are_compatible("C", "G")
    assert not are_compatible("A", "A")
    assert not are_compatible("G", "A")


def test_degenerate_codes_are_set_compatible():
    # R=A/G，其互补集合为 T/C = Y
    assert are_compatible("R", "Y")
    assert are_compatible("Y", "R")
    # R 与 W(A/T) 在 A 上有交集
    assert are_compatible("R", "W")
    # S=G/C 与 W=A/T 无交集
    assert not are_compatible("S", "W")
    # N 与任何碱基都可能互补
    for code in "ATCGRYSWKMBDHVN":
        assert are_compatible("N", code)
        assert are_compatible(code, "N")
    # B=C/G/T 互补集合为 G/C/A = V
    assert are_compatible("B", "V")


def test_validate_iupac():
    assert validate_iupac(" acgtn ") == "ACGTN"
    assert validate_iupac("rYswkmbdhv") == "RYSWKMBDHV"
    for bad in ["", "ACXG", "AC-GT", "AC GT"]:
        with pytest.raises(ValueError):
            validate_iupac(bad)


# ---------- 对齐枚举 ----------

def test_all_shifts_enumerated():
    r = analyze_pair("P", "AAAAA", "Q", "CCCGGG", 4, 3)
    # a + b - 1 = 5 + 6 - 1 = 10
    assert r.alignments_checked == 10
    aligns = [
        _enumerate_alignment("AAAAA", "CCCGGG", d)
        for d in range(-5, 5)
    ]
    assert [a.overlap_len for a in aligns] == [1, 2, 3, 4, 5, 5, 4, 3, 2, 1]


def test_full_complement_flush():
    r = analyze_pair("P1", "AAAAA", "P2", "TTTTT", 4, 4)
    assert r.best_shift == 0
    assert r.max_complementary_run == 5
    assert r.top3_run == 5
    assert r.bottom3_run == 5
    assert r.terminal_run == 5
    assert r.risk_level == "high"
    assert r.overlap_len == 5


def test_three_prime_tip_dimer():
    # 两条引物各自 3′端 5 个碱基互补
    r = analyze_pair(
        "F", "GGGGGAAAAA", "R", "GGGGGTTTTT", 4, 3
    )
    assert r.max_complementary_run == 5
    assert r.top3_run == 5
    assert r.bottom3_run == 5
    assert r.risk_level == "high"
    assert "延伸风险" in r.risk_reasons[0]


def test_extreme_shifts_geometry():
    left = _enumerate_alignment("AAAAA", "TTTTT", -4)
    right = _enumerate_alignment("AAAAA", "TTTTT", 4)
    assert left.overlap_len == right.overlap_len == 1
    assert left.top3_run == left.bottom3_run == 0
    assert right.top3_run == right.bottom3_run == 1
    # 渲染行等宽
    assert len(left.top_line) == len(left.bottom_line) == len(left.match_line)
    assert len(right.top_line) == len(right.bottom_line) == len(right.match_line)
    # 右极端错位时唯一配对列显示 T4 | T(b-1)
    assert left.match_line.strip() == "|"
    assert right.match_line.strip() == "|"


def test_medium_risk_stable_run_without_tip():
    # 内部 4 连续互补，3′端配对均错配 -> medium 而非 high
    r = analyze_pair(
        "P", "GGCCCCGG", "Q", "TTGGGGTT", 4, 4
    )
    assert r.max_complementary_run == 4
    assert r.terminal_run < 4
    assert r.risk_level == "medium"


def test_low_risk():
    r = analyze_pair("P", "AAAAAAAA", "Q", "AAAAAAAA", 4, 3)
    assert r.max_complementary_run == 0
    assert r.risk_level == "low"


def test_self_dimer_flag_and_count():
    r = analyze_pair("P", "AAAATTTT", "P", "AAAATTTT", 3, 3, is_self_dimer=True)
    assert r.is_self_dimer
    assert r.alignments_checked == 15  # 2*8-1
    assert r.name_top == r.name_bottom


def test_tie_break_is_unique_and_stable():
    # 自二聚体完全不互补：所有对齐风险键相同，固定偏移规则应选 |d| 最小的 d=0
    r = analyze_pair("P", "AAAAAAAA", "P", "AAAAAAAA", 2, 2, is_self_dimer=True)
    assert r.best_shift == 0

    # 多次计算结果一致
    again = analyze_pair(
        "P", "AAAAAAAA", "P", "AAAAAAAA", 2, 2, is_self_dimer=True
    )
    assert r.best_shift == again.best_shift == 0


def test_tie_break_key_smallest_abs_shift():
    # 全部 N 的自二聚体：每个错位的风险键完全相同
    # （每个重叠列都互补、tip run 等于重叠长度），固定偏移规则应选 |d| 最小的 d=0
    r = analyze_pair("P", "NNN", "P", "NNN", 1, 1, is_self_dimer=True)
    assert r.best_shift == 0
    assert r.max_complementary_run == 3  # d=0 时重叠 3 列全互补


def test_tie_break_synthetic_positive_preferred():
    # 直接验证同分键：风险键相同时正偏移优先于负偏移，小 |shift| 优先
    from app.dimer import Alignment, _tie_key

    def aln(d: int) -> Alignment:
        return Alignment(
            shift=d, overlap_len=2, max_complementary_run=1,
            top3_run=0, bottom3_run=0,
            top_line="", match_line="", bottom_line="",
        )

    assert min(aln(1), aln(-1), key=_tie_key).shift == 1
    assert min(aln(0), aln(2), key=_tie_key).shift == 0
    assert min(aln(2), aln(-2), key=_tie_key).shift == 2  # 正优先
    assert min(aln(1), aln(0), key=_tie_key).shift == 0   # |d| 优先于符号
    # 全序：任意键可排序且稳定
    ordered = sorted([aln(-2), aln(2), aln(1), aln(-1), aln(0)], key=_tie_key)
    assert [a.shift for a in ordered] == [0, 1, -1, 2, -2]
