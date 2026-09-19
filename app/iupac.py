"""IUPAC 简并碱基集合与互补判断。

每个 IUPAC 简并码映射为其代表的碱基集合；两条链上的碱基“可能互补”
当且仅当其中一个集合与另一个集合的 Watson-Crick 互补集合存在交集。
"""

from __future__ import annotations

from functools import lru_cache

# IUPAC nucleotide ambiguity codes（仅非缺口碱基；- / . 等缺口不接受）
IUPAC_BASE_SETS: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "T": frozenset("T"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "R": frozenset("AG"),  # puRine
    "Y": frozenset("CT"),  # pYrimidine
    "S": frozenset("GC"),  # Strong
    "W": frozenset("AT"),  # Weak
    "K": frozenset("GT"),  # Keto
    "M": frozenset("AC"),  # aMino
    "B": frozenset("CGT"),  # not A
    "D": frozenset("AGT"),  # not C
    "H": frozenset("ACT"),  # not G
    "V": frozenset("ACG"),  # not T
    "N": frozenset("ACGT"),  # aNy
}

_VALID_CHARS = frozenset(IUPAC_BASE_SETS)

# Watson-Crick 互补
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def normalize_sequence(sequence: str) -> str:
    """去除首尾空白并转为大写；不做内部空白或非法字符的容忍处理。"""
    return sequence.strip().upper()


def validate_iupac(sequence: str) -> str:
    """校验是否全部由 IUPAC 简并碱基组成，返回规范化后的大写序列。

    非法字符抛出 ValueError，由 API 层转换为 422。
    """
    seq = normalize_sequence(sequence)
    if not seq:
        raise ValueError("引物序列不能为空")
    bad = sorted({ch for ch in seq if ch not in _VALID_CHARS})
    if bad:
        raise ValueError(f"含非法 IUPAC 碱基: {''.join(bad)}")
    return seq


@lru_cache(maxsize=1)
def _compatibility_matrix() -> dict[tuple[str, str], bool]:
    """预计算 15x15 简并码“可能互补”矩阵。"""
    matrix: dict[tuple[str, str], bool] = {}
    for code_a, set_a in IUPAC_BASE_SETS.items():
        complement_a = {_COMPLEMENT[b] for b in set_a}
        for code_b, set_b in IUPAC_BASE_SETS.items():
            matrix[(code_a, code_b)] = bool(complement_a & set_b)
    return matrix


def are_compatible(base_top: str, base_bottom: str) -> bool:
    """两个 IUPAC 碱基是否可能 Watson-Crick 互补（集合交集非空）。"""
    return _compatibility_matrix()[(base_top, base_bottom)]
