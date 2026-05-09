"""
residue_props.py

本文件存放标准氨基酸的基础理化性质。

为什么单独写这个文件？
因为残基理化性质会在多个地方用到，例如：
1. 提取 residue-level features；
2. 聚合 nanopore-level features；
3. 后续构建残基图；
4. 后续筛选关键孔道残基。

当前第一版只使用比较基础、稳定的性质：
- 电荷 charge；
- Kyte-Doolittle 疏水性 hydrophobicity；
- 是否芳香族 is_aromatic；
- 是否极性 is_polar。

注意：
这里的电荷是简化近似值，不是严格的 pH 条件下质子化状态。
例如 HIS 在不同 pH 下可能带正电，也可能近似中性。
第一版为了工程可运行，先给 HIS 一个弱正电近似值 0.1。
后续如果需要结合 pH，可以单独建立 pH-adjusted charge 特征。
"""

from typing import Dict, Optional


# 20 种标准氨基酸三字母缩写
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}


# 简化残基电荷表
# ASP / GLU 通常带负电
# LYS / ARG 通常带正电
# HIS 受 pH 影响较大，这里给一个弱正电近似
RESIDUE_CHARGE: Dict[str, float] = {
    "ASP": -1.0,
    "GLU": -1.0,
    "LYS": 1.0,
    "ARG": 1.0,
    "HIS": 0.1,
}


# Kyte-Doolittle hydrophobicity scale
# 数值越大，越疏水；数值越小，越亲水
HYDROPHOBICITY: Dict[str, float] = {
    "ILE": 4.5,
    "VAL": 4.2,
    "LEU": 3.8,
    "PHE": 2.8,
    "CYS": 2.5,
    "MET": 1.9,
    "ALA": 1.8,
    "GLY": -0.4,
    "THR": -0.7,
    "SER": -0.8,
    "TRP": -0.9,
    "TYR": -1.3,
    "PRO": -1.6,
    "HIS": -3.2,
    "GLU": -3.5,
    "GLN": -3.5,
    "ASP": -3.5,
    "ASN": -3.5,
    "LYS": -3.9,
    "ARG": -4.5,
}


# 芳香族氨基酸
# HIS 是否算芳香族取决于定义，这里把 HIS 也纳入芳香杂环类
AROMATIC_RESIDUES = {
    "PHE",
    "TYR",
    "TRP",
    "HIS",
}


# 极性或带电残基
POLAR_RESIDUES = {
    "SER",
    "THR",
    "ASN",
    "GLN",
    "TYR",
    "CYS",
    "LYS",
    "ARG",
    "HIS",
    "ASP",
    "GLU",
}


def normalize_residue_name(resname: str) -> str:
    """
    标准化残基名称。

    参数
    ----
    resname : str
        Biopython 读取到的残基名称，例如 "ALA"。

    返回
    ----
    str
        去除空格并转换为大写后的残基名称。
    """

    return resname.strip().upper()


def is_standard_amino_acid(resname: str) -> bool:
    """
    判断残基名称是否属于 20 种标准氨基酸。

    参数
    ----
    resname : str
        三字母残基名称。

    返回
    ----
    bool
        True 表示标准氨基酸；
        False 表示非标准残基、水分子、配体等。
    """

    resname = normalize_residue_name(resname)
    return resname in STANDARD_AMINO_ACIDS


def get_residue_charge(resname: str) -> float:
    """
    获取残基的简化电荷。

    参数
    ----
    resname : str
        三字母残基名称。

    返回
    ----
    float
        简化电荷值。
        如果不在 RESIDUE_CHARGE 表中，默认为 0.0。
    """

    resname = normalize_residue_name(resname)
    return RESIDUE_CHARGE.get(resname, 0.0)


def get_hydrophobicity(resname: str) -> Optional[float]:
    """
    获取残基的 Kyte-Doolittle 疏水性数值。

    参数
    ----
    resname : str
        三字母残基名称。

    返回
    ----
    float | None
        疏水性数值。
        如果不是标准氨基酸，则返回 None。
    """

    resname = normalize_residue_name(resname)
    return HYDROPHOBICITY.get(resname, None)


def is_aromatic(resname: str) -> int:
    """
    判断是否为芳香族残基。

    返回 int 而不是 bool，是为了后续更方便保存到 CSV：
    1 表示是；
    0 表示否。
    """

    resname = normalize_residue_name(resname)
    return int(resname in AROMATIC_RESIDUES)


def is_polar(resname: str) -> int:
    """
    判断是否为极性或带电残基。

    返回
    ----
    int
        1 表示极性/带电；
        0 表示非极性。
    """

    resname = normalize_residue_name(resname)
    return int(resname in POLAR_RESIDUES)