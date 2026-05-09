"""
pore.py

本文件负责与纳米孔孔道相关的处理。

当前阶段主要实现：
1. 读取 residue_features.csv；
2. 基于几何规则添加 radial_distance 等特征；
3. 根据径向距离阈值筛选候选孔道内壁残基；
4. 保存 inner_candidate_residues.csv。

重要说明
----
当前方法是第一版几何近似方法：

    radial_distance <= inner_radius_threshold

它不能严格等同于专业孔道识别工具给出的 pore-lining residues。

但它的优点是：
1. 简单；
2. 可控；
3. 适合先跑通工程流程；
4. 后续可以替换为 HOLE / CHAP / MOLEonline 的结果。
"""

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from npstructfeat.geometry import add_basic_geometry_features
from npstructfeat.io import get_output_file, resolve_path


def get_residue_feature_file(config: Dict[str, Any]) -> Path:
    """
    根据配置文件获取 residue_features.csv 的路径。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    residue_feature_file : Path
        残基级特征表路径。

    文件命名规则
    ----
    如果 pdb_id = 5JZT，则默认读取：
        data/processed/residue_features/5JZT_residue_features.csv
    """

    residue_feature_file = get_output_file(
        config=config,
        output_dir_key="residue_feature_dir",
        suffix="residue_features.csv",
    )

    if not residue_feature_file.exists():
        raise FileNotFoundError(
            f"残基级特征表不存在: {residue_feature_file}\n"
            "请先运行 scripts/02_extract_residue_features.py"
        )

    return residue_feature_file


def load_residue_features(config: Dict[str, Any]) -> pd.DataFrame:
    """
    读取 residue_features.csv。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    df : pandas.DataFrame
        残基级特征表。
    """

    residue_feature_file = get_residue_feature_file(config)

    df = pd.read_csv(residue_feature_file)

    return df


def select_inner_candidate_residues(config: Dict[str, Any]) -> pd.DataFrame:
    """
    筛选候选孔道内壁残基。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    inner_df : pandas.DataFrame
        候选孔道内壁残基表。

    使用的配置字段
    ----
    config["pore"]["axis_mode"]
        当前只支持 "z_axis"。

    config["pore"]["center_mode"]
        当前支持 "xy_mean" 或 "xy_median"。

    config["pore"]["inner_radius_threshold"]
        径向距离阈值，单位 Å。

    当前筛选规则
    ----
    1. 假设孔道轴方向是 z 轴；
    2. 在 x-y 平面估计中心点；
    3. 计算每个残基到中心轴的径向距离；
    4. 保留 radial_distance <= threshold 的残基。
    """

    residue_df = load_residue_features(config)

    pore_config = config.get("pore", {})

    axis_mode = pore_config.get("axis_mode", "z_axis")
    center_mode = pore_config.get("center_mode", "xy_mean")
    threshold = float(pore_config.get("inner_radius_threshold", 20.0))

    if axis_mode != "z_axis":
        raise ValueError(
            f"当前第一版仅支持 axis_mode='z_axis'，但得到: {axis_mode}\n"
            "如果结构未沿 z 轴对齐，需要后续增加 PCA 主轴对齐或专业孔道识别工具。"
        )

    # 添加 radial_distance、theta、z_norm 等几何特征
    residue_geo_df = add_basic_geometry_features(
        residue_df=residue_df,
        center_mode=center_mode,
    )

    # 按径向距离筛选候选孔道内壁残基
    inner_df = residue_geo_df[
        residue_geo_df["radial_distance"] <= threshold
    ].copy()

    # 标记为候选内壁残基
    inner_df["is_inner_candidate"] = 1
    inner_df["inner_radius_threshold"] = threshold
    inner_df["axis_mode"] = axis_mode
    inner_df["center_mode"] = center_mode

    # 为了阅读方便，按 z 坐标和链 ID 排序
    inner_df = inner_df.sort_values(
        by=["z", "chain_id", "residue_number"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    return inner_df


def save_inner_candidate_residues(config: Dict[str, Any]) -> Path:
    """
    筛选并保存候选孔道内壁残基表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_file : Path
        保存后的 CSV 文件路径。
    """

    inner_df = select_inner_candidate_residues(config)

    output_file = get_output_file(
        config=config,
        output_dir_key="inner_residue_dir",
        suffix="inner_candidate_residues.csv",
    )

    inner_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file


def summarize_inner_candidates(inner_df: pd.DataFrame) -> pd.DataFrame:
    """
    对候选孔道内壁残基进行链级别统计。

    参数
    ----
    inner_df : pandas.DataFrame
        候选孔道内壁残基表。

    返回
    ----
    summary_df : pandas.DataFrame
        每条链的候选内壁残基数量和理化性质统计。

    统计内容
    ----
    - inner_candidate_count
    - mean_radial_distance
    - min_radial_distance
    - max_radial_distance
    - mean_z
    - min_z
    - max_z
    - net_charge
    - mean_hydrophobicity
    - aromatic_count
    - polar_count
    """

    if len(inner_df) == 0:
        return pd.DataFrame()

    summary_df = (
        inner_df.groupby("chain_id")
        .agg(
            inner_candidate_count=("residue_number", "count"),
            mean_radial_distance=("radial_distance", "mean"),
            min_radial_distance=("radial_distance", "min"),
            max_radial_distance=("radial_distance", "max"),
            mean_z=("z", "mean"),
            min_z=("z", "min"),
            max_z=("z", "max"),
            net_charge=("charge", "sum"),
            mean_hydrophobicity=("hydrophobicity", "mean"),
            aromatic_count=("is_aromatic", "sum"),
            polar_count=("is_polar", "sum"),
        )
        .reset_index()
    )

    return summary_df