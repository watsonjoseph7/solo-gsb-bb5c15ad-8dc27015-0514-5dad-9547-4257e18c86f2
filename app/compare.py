"""批次间配对风险对比：稳定配对键、风险升降分类。

配对键
------
按引物**名称**构造稳定配对键，与批次内提交顺序无关：

* 自二聚体：保留单名，键为 ``(name, name)``；
* 异二聚体：忽略批次内提交顺序，两名按字典序排列，键为 ``(小名, 大名)``。

同名引物在两批次间序列发生变化时，仍按同一配对键比较，并在明细中标记
``sequence_changed``。

风险等级全序为 ``low < medium < high``。两批次共有的配对按键比较等级：
上升为 escalated、下降为 deescalated、相同不计入明细；仅存在于 target
的配对为 introduced，仅存在于 base 的为 resolved。四类明细各自按配对键
（即配对名字典序）确定性排序，相同输入必得相同输出。
"""

from __future__ import annotations

from dataclasses import dataclass

# 风险等级全序：low < medium < high
RISK_ORDER = {"low": 0, "medium": 1, "high": 2}


@dataclass(frozen=True)
class PairDiff:
    """一个配对键在两批次间的风险变化明细。"""

    name_a: str                     # 配对键第一名称（异二聚体取字典序较小者）
    name_b: str                     # 配对键第二名称（自二聚体与 name_a 相同）
    is_self_dimer: bool
    sequence_changed: bool          # 配对中任一同名引物在两批次间序列不同
    base_risk_level: str | None     # base 批次风险等级（introduced 时为 None）
    target_risk_level: str | None   # target 批次风险等级（resolved 时为 None）

    @property
    def pair_key(self) -> str:
        """配对键的可读形式：自二聚体为单名，异二聚体为 ``A|B``。"""
        if self.is_self_dimer:
            return self.name_a
        return f"{self.name_a}|{self.name_b}"


@dataclass(frozen=True)
class CompareResult:
    """两批次对比结果：四类明细（各按配对键升序）与共有配对计数。"""

    escalated: tuple[PairDiff, ...]
    deescalated: tuple[PairDiff, ...]
    introduced: tuple[PairDiff, ...]
    resolved: tuple[PairDiff, ...]
    common_pair_count: int          # 两批次共有的配对键数
    unchanged_pair_count: int       # 共有配对中风险等级不变的数量


def pair_key(name_top: str, name_bottom: str, is_self_dimer: bool) -> tuple[str, str]:
    """构造稳定配对键（见模块 docstring）。"""
    if is_self_dimer or name_top == name_bottom:
        return (name_top, name_top)
    if name_top < name_bottom:
        return (name_top, name_bottom)
    return (name_bottom, name_top)


def compare_pair_risks(
    *,
    base_pairs: list[tuple[str, str, bool, str]],
    target_pairs: list[tuple[str, str, bool, str]],
    base_sequences: dict[str, str],
    target_sequences: dict[str, str],
) -> CompareResult:
    """按配对键对比两批次的风险等级。

    base_pairs / target_pairs 元素为
    ``(name_top, name_bottom, is_self_dimer, risk_level)``；
    base_sequences / target_sequences 为 ``引物名 -> 序列``。
    """
    base_map = {
        pair_key(top, bottom, is_self): level
        for top, bottom, is_self, level in base_pairs
    }
    target_map = {
        pair_key(top, bottom, is_self): level
        for top, bottom, is_self, level in target_pairs
    }

    # 两批次都存在但序列不同的引物名
    changed_names = {
        name
        for name, seq in base_sequences.items()
        if name in target_sequences and target_sequences[name] != seq
    }

    escalated: list[PairDiff] = []
    deescalated: list[PairDiff] = []
    introduced: list[PairDiff] = []
    resolved: list[PairDiff] = []
    common = unchanged = 0

    # 按配对键（名称元组字典序）遍历，四类明细天然确定性排序
    for key in sorted(set(base_map) | set(target_map)):
        name_a, name_b = key
        base_level = base_map.get(key)
        target_level = target_map.get(key)
        diff = PairDiff(
            name_a=name_a,
            name_b=name_b,
            is_self_dimer=name_a == name_b,
            sequence_changed=name_a in changed_names or name_b in changed_names,
            base_risk_level=base_level,
            target_risk_level=target_level,
        )
        if base_level is None:
            introduced.append(diff)
        elif target_level is None:
            resolved.append(diff)
        else:
            common += 1
            if RISK_ORDER[target_level] > RISK_ORDER[base_level]:
                escalated.append(diff)
            elif RISK_ORDER[target_level] < RISK_ORDER[base_level]:
                deescalated.append(diff)
            else:
                unchanged += 1

    return CompareResult(
        escalated=tuple(escalated),
        deescalated=tuple(deescalated),
        introduced=tuple(introduced),
        resolved=tuple(resolved),
        common_pair_count=common,
        unchanged_pair_count=unchanged,
    )
