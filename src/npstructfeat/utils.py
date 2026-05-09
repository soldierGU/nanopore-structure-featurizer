"""
utils.py

本文件存放一些通用工具函数。

这些函数不专属于结构解析，也不专属于残基特征提取，
而是整个工程都可能会用到的辅助函数。

例如：
1. 打印配置摘要；
2. 标准化 PDB ID；
3. 检查配置字段是否存在。
"""

from typing import Any, Dict


def normalize_pdb_id(pdb_id: str) -> str:
    """
    标准化 PDB ID。

    PDB ID 通常使用大写，例如：
    5JZT
    9FM6
    1PRE

    参数
    ----
    pdb_id : str
        原始 PDB ID。

    返回
    ----
    normalized_id : str
        大写后的 PDB ID。
    """

    return pdb_id.strip().upper()


def print_config_summary(config: Dict[str, Any]) -> None:
    """
    打印当前配置摘要。

    这个函数主要用于初学阶段检查配置是否读取正确。

    参数
    ----
    config : dict
        从 YAML 读取到的配置字典。
    """

    input_config = config.get("input", {})
    pore_config = config.get("pore", {})
    structure_config = config.get("structure", {})
    output_config = config.get("output", {})

    print("=" * 60)
    print("Nanopore Structure Featurizer 配置摘要")
    print("=" * 60)

    print("[Input]")
    print(f"PDB ID:          {input_config.get('pdb_id')}")
    print(f"Nanopore ID:     {input_config.get('nanopore_id')}")
    print(f"Structure file:  {input_config.get('structure_file')}")
    print(f"File format:     {input_config.get('file_format')}")
    print(f"Assembly type:   {input_config.get('assembly_type')}")

    print("\n[Structure]")
    print(f"Model index:     {structure_config.get('use_model_index')}")
    print(f"Standard only:   {structure_config.get('use_standard_residues_only')}")
    print(f"Use CA only:     {structure_config.get('use_ca_only')}")

    print("\n[Pore]")
    print(f"Axis mode:       {pore_config.get('axis_mode')}")
    print(f"Center mode:     {pore_config.get('center_mode')}")
    print(f"Inner threshold: {pore_config.get('inner_radius_threshold')} Å")

    print("\n[Output]")
    for key, value in output_config.items():
        print(f"{key}: {value}")

    print("=" * 60)


def require_config_keys(config: Dict[str, Any], required_keys: list[str]) -> None:
    """
    检查配置文件是否包含必要的一级字段。

    参数
    ----
    config : dict
        YAML 配置字典。

    required_keys : list[str]
        必须存在的一级字段。
        例如：
        ["input", "structure", "pore", "output"]

    如果缺少字段，会直接抛出 KeyError。

    为什么需要这个函数？
    因为配置文件一旦写错，后面的代码可能会报出很难理解的错误。
    提前检查可以让错误更清楚。
    """

    for key in required_keys:
        if key not in config:
            raise KeyError(f"配置文件缺少必要字段: {key}")