"""
features.py

本文件负责从 PDB/mmCIF 结构中提取特征。

当前阶段主要实现：
1. 提取 residue-level features；
2. 每个标准氨基酸残基生成一行；
3. 使用 Cα 原子坐标代表残基空间位置；
4. 加入残基基础理化性质；
5. 保存为 CSV 文件。

后续还会在这个文件中继续加入：
- nanopore-level 聚合特征；
- 孔道内壁残基聚合特征；
- 结构统计特征等。
"""

from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio.PDB.Residue import Residue

from npstructfeat.io import get_output_file
from npstructfeat.parser import get_selected_model, is_standard_residue, load_structure
from npstructfeat.residue_props import (
    get_hydrophobicity,
    get_residue_charge,
    is_aromatic,
    is_polar,
    normalize_residue_name,
)


def get_residue_number(residue: Residue) -> int:
    """
    获取残基编号。

    Biopython 中 residue.id 通常是一个三元组：
        residue.id = (hetero_flag, sequence_identifier, insertion_code)

    例如普通残基：
        (" ", 238, " ")

    其中 residue.id[1] 就是残基编号。

    参数
    ----
    residue : Bio.PDB.Residue.Residue
        Biopython 残基对象。

    返回
    ----
    int
        残基编号。
    """

    return residue.id[1]


def get_insertion_code(residue: Residue) -> str:
    """
    获取残基插入码 insertion code。

    大多数结构中 insertion code 是空格。
    有些 PDB 为了处理编号插入，会使用 A、B 等插入码。

    例如：
        100
        100A
        100B

    第一版先记录下来，方便后续排查。
    """

    insertion_code = residue.id[2]

    if insertion_code == " ":
        return ""

    return str(insertion_code).strip()


def extract_ca_coordinates(residue: Residue):
    """
    提取残基的 Cα 坐标。

    参数
    ----
    residue : Bio.PDB.Residue.Residue
        残基对象。

    返回
    ----
    tuple[float, float, float] | None
        如果残基有 CA 原子，返回 (x, y, z)；
        如果没有 CA 原子，返回 None。

    为什么用 CA？
    对于残基级别建模，Cα 坐标常被用来表示残基的空间位置。
    这比使用所有原子更简单，也更适合第一版工程。
    """

    if "CA" not in residue:
        return None

    ca_atom = residue["CA"]
    x, y, z = ca_atom.coord

    return float(x), float(y), float(z)


def extract_residue_features(config: Dict[str, Any]) -> pd.DataFrame:
    """
    从结构文件中提取标准氨基酸残基特征。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    df : pandas.DataFrame
        残基级别特征表。

    表格每一行对应一个标准氨基酸残基。

    当前提取的特征包括：
    - pdb_id
    - nanopore_id
    - model_index
    - chain_id
    - residue_number
    - insertion_code
    - residue_name
    - x, y, z
    - charge
    - hydrophobicity
    - is_aromatic
    - is_polar
    - has_ca

    注意：
    默认只保留标准氨基酸残基。
    如果某个标准残基缺失 CA 原子，则仍然记录该残基，
    但 x/y/z 会记为 None，has_ca = 0。
    后续如果 use_ca_only = true，可以再过滤掉这些残基。
    """

    structure = load_structure(config)

    model_index = config.get("structure", {}).get("use_model_index", 0)
    model = get_selected_model(structure, model_index=model_index)

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")

    use_standard_only = config.get("structure", {}).get(
        "use_standard_residues_only", True
    )
    use_ca_only = config.get("structure", {}).get("use_ca_only", True)

    rows: List[Dict[str, Any]] = []

    for chain in model:
        chain_id = chain.id

        for residue in chain:
            # 第一版默认只处理标准氨基酸残基
            if use_standard_only and not is_standard_residue(residue):
                continue

            resname = normalize_residue_name(residue.resname)

            # 提取 Cα 坐标
            ca_coord = extract_ca_coordinates(residue)

            if ca_coord is None:
                x, y, z = None, None, None
                has_ca = 0
            else:
                x, y, z = ca_coord
                has_ca = 1

            # 如果配置要求只使用有 CA 的残基，则跳过缺失 CA 的残基
            if use_ca_only and has_ca == 0:
                continue

            row = {
                "pdb_id": pdb_id,
                "nanopore_id": nanopore_id,
                "model_index": model_index,
                "chain_id": chain_id,
                "residue_number": get_residue_number(residue),
                "insertion_code": get_insertion_code(residue),
                "residue_name": resname,
                "x": x,
                "y": y,
                "z": z,
                "has_ca": has_ca,
                "charge": get_residue_charge(resname),
                "hydrophobicity": get_hydrophobicity(resname),
                "is_aromatic": is_aromatic(resname),
                "is_polar": is_polar(resname),
            }

            rows.append(row)

    df = pd.DataFrame(rows)

    # 统一列顺序，便于查看和后续读取
    columns = [
        "pdb_id",
        "nanopore_id",
        "model_index",
        "chain_id",
        "residue_number",
        "insertion_code",
        "residue_name",
        "x",
        "y",
        "z",
        "has_ca",
        "charge",
        "hydrophobicity",
        "is_aromatic",
        "is_polar",
    ]

    df = df[columns]

    return df


def save_residue_features(config: Dict[str, Any]) -> Path:
    """
    提取并保存 residue-level features。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_file : Path
        保存后的 CSV 文件路径。
    """

    df = extract_residue_features(config)

    output_file = get_output_file(
        config=config,
        output_dir_key="residue_feature_dir",
        suffix="residue_features.csv",
    )

    df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file

def get_inner_candidate_file(config: Dict[str, Any]) -> Path:
    """
    根据配置文件获取候选孔道内壁残基表路径。

    默认路径规则：
        data/processed/inner_residues/<PDB_ID>_inner_candidate_residues.csv

    例如：
        data/processed/inner_residues/5JZT_inner_candidate_residues.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    inner_file : Path
        候选孔道内壁残基 CSV 文件路径。

    注意
    ----
    运行本函数前，需要先执行：
        scripts/03_select_inner_residues.py
    """

    inner_file = get_output_file(
        config=config,
        output_dir_key="inner_residue_dir",
        suffix="inner_candidate_residues.csv",
    )

    if not inner_file.exists():
        raise FileNotFoundError(
            f"候选孔道内壁残基表不存在: {inner_file}\n"
            "请先运行 scripts/03_select_inner_residues.py"
        )

    return inner_file


def build_nanopore_structure_features(config: Dict[str, Any]) -> pd.DataFrame:
    """
    将候选孔道内壁残基表聚合为 nanopore-level 结构特征表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    feature_df : pandas.DataFrame
        纳米孔整体结构特征表。

    表格含义
    ----
    每一行对应一个纳米孔结构。

    当前你的数据中只有：
        5JZT / aerolysin_WT

    所以输出表只有一行。

    主要特征包括：
    1. 候选孔道内壁残基数量；
    2. 内壁残基净电荷；
    3. 内壁平均疏水性；
    4. 芳香族残基数量和比例；
    5. 极性残基数量和比例；
    6. z 方向覆盖范围；
    7. 径向距离统计；
    8. 每条链候选残基数量的均值和标准差。

    注意
    ----
    这些特征来自几何粗筛结果。
    因此更准确的命名是：
        candidate inner residues features
    而不是严格的 confirmed pore-lining residue features。
    """

    inner_file = get_inner_candidate_file(config)
    inner_df = pd.read_csv(inner_file)

    if len(inner_df) == 0:
        raise ValueError(
            "候选孔道内壁残基表为空，无法构建 nanopore-level features。"
        )

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")
    structure_file = config["input"].get("structure_file", "")
    assembly_type = config["input"].get("assembly_type", "")
    file_format = config["input"].get("file_format", "")

    threshold = float(config.get("pore", {}).get("inner_radius_threshold", 20.0))
    axis_mode = config.get("pore", {}).get("axis_mode", "z_axis")
    center_mode = config.get("pore", {}).get("center_mode", "xy_mean")

    # 总残基数：来自 residue_features 表
    # 用于计算候选内壁残基占比
    residue_feature_file = get_output_file(
        config=config,
        output_dir_key="residue_feature_dir",
        suffix="residue_features.csv",
    )

    if residue_feature_file.exists():
        residue_df = pd.read_csv(residue_feature_file)
        total_residue_count = len(residue_df)
    else:
        # 理论上不应该发生；保留兜底逻辑
        total_residue_count = None

    inner_residue_count = len(inner_df)

    if total_residue_count is not None and total_residue_count > 0:
        inner_residue_ratio = inner_residue_count / total_residue_count
    else:
        inner_residue_ratio = None

    # 按链统计候选残基数量，用于检查七聚体对称性
    chain_counts = inner_df.groupby("chain_id")["residue_number"].count()

    # z 方向范围：粗略表示候选孔道区域覆盖长度
    z_min = float(inner_df["z"].min())
    z_max = float(inner_df["z"].max())
    z_range = z_max - z_min

    # 径向距离统计
    radial_min = float(inner_df["radial_distance"].min())
    radial_max = float(inner_df["radial_distance"].max())
    radial_mean = float(inner_df["radial_distance"].mean())
    radial_std = float(inner_df["radial_distance"].std())

    # 理化性质聚合
    inner_net_charge = float(inner_df["charge"].sum())
    inner_mean_charge = float(inner_df["charge"].mean())

    inner_mean_hydrophobicity = float(inner_df["hydrophobicity"].mean())
    inner_std_hydrophobicity = float(inner_df["hydrophobicity"].std())

    aromatic_count = int(inner_df["is_aromatic"].sum())
    aromatic_ratio = float(inner_df["is_aromatic"].mean())

    polar_count = int(inner_df["is_polar"].sum())
    polar_ratio = float(inner_df["is_polar"].mean())

    # 正负电残基数量
    positive_count = int((inner_df["charge"] > 0).sum())
    negative_count = int((inner_df["charge"] < 0).sum())
    neutral_count = int((inner_df["charge"] == 0).sum())

    # 残基种类数
    residue_type_count = int(inner_df["residue_name"].nunique())

    # 中心轴估计值
    # 这些列由 03_select_inner_residues.py 添加
    center_x = float(inner_df["center_x"].iloc[0])
    center_y = float(inner_df["center_y"].iloc[0])

    features = {
        # 基础信息
        "nanopore_id": nanopore_id,
        "pdb_id": pdb_id,
        "structure_file": structure_file,
        "file_format": file_format,
        "assembly_type": assembly_type,

        # 方法参数
        "axis_mode": axis_mode,
        "center_mode": center_mode,
        "inner_radius_threshold_A": threshold,
        "center_x": center_x,
        "center_y": center_y,

        # 残基数量信息
        "total_residue_count": total_residue_count,
        "inner_candidate_residue_count": inner_residue_count,
        "inner_candidate_residue_ratio": inner_residue_ratio,
        "chain_count": int(inner_df["chain_id"].nunique()),
        "inner_candidate_count_per_chain_mean": float(chain_counts.mean()),
        "inner_candidate_count_per_chain_std": float(chain_counts.std()),
        "inner_candidate_count_per_chain_min": int(chain_counts.min()),
        "inner_candidate_count_per_chain_max": int(chain_counts.max()),

        # z 方向孔道区域范围
        "z_min_A": z_min,
        "z_max_A": z_max,
        "pore_region_length_approx_A": z_range,

        # 径向距离统计
        "radial_distance_min_A": radial_min,
        "radial_distance_max_A": radial_max,
        "radial_distance_mean_A": radial_mean,
        "radial_distance_std_A": radial_std,

        # 电荷相关特征
        "inner_net_charge": inner_net_charge,
        "inner_mean_charge": inner_mean_charge,
        "inner_positive_residue_count": positive_count,
        "inner_negative_residue_count": negative_count,
        "inner_neutral_residue_count": neutral_count,

        # 疏水性相关特征
        "inner_mean_hydrophobicity": inner_mean_hydrophobicity,
        "inner_std_hydrophobicity": inner_std_hydrophobicity,

        # 芳香族/极性残基特征
        "inner_aromatic_count": aromatic_count,
        "inner_aromatic_ratio": aromatic_ratio,
        "inner_polar_count": polar_count,
        "inner_polar_ratio": polar_ratio,

        # 残基组成复杂度
        "inner_residue_type_count": residue_type_count,
    }

    feature_df = pd.DataFrame([features])

    return feature_df


def save_nanopore_structure_features(config: Dict[str, Any]) -> Path:
    """
    构建并保存 nanopore-level 结构特征表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_file : Path
        保存后的 CSV 文件路径。
    """

    feature_df = build_nanopore_structure_features(config)

    output_file = get_output_file(
        config=config,
        output_dir_key="nanopore_feature_dir",
        suffix="nanopore_structure_features.csv",
    )

    feature_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file