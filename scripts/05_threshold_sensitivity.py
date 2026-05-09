"""
05_threshold_sensitivity.py

本脚本用于分析 inner_radius_threshold 对候选孔道内壁残基筛选结果的影响。

为什么需要这个脚本？
当前我们使用的候选孔道内壁残基筛选规则是：

    radial_distance <= inner_radius_threshold

其中 inner_radius_threshold 是人为设定的几何阈值，例如 20 Å。
如果不同阈值下提取到的结构特征变化很剧烈，说明该方法对阈值敏感；
如果特征变化比较平滑，说明该几何粗筛方法作为第一版工程特征较稳定。

运行方式：
    python scripts/05_threshold_sensitivity.py

或者指定配置文件：
    python scripts/05_threshold_sensitivity.py --config config/default.yaml

也可以指定阈值：
    python scripts/05_threshold_sensitivity.py --thresholds 15 18 20 22 25

输入文件：
    data/processed/residue_features/5JZT_residue_features.csv

输出文件：
    data/processed/nanopore_features/5JZT_threshold_sensitivity.csv

注意：
本脚本不会覆盖：
    5JZT_inner_candidate_residues.csv
    5JZT_nanopore_structure_features.csv

它只额外输出 threshold sensitivity 分析表。
"""

import argparse
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from npstructfeat.geometry import add_basic_geometry_features
from npstructfeat.io import get_output_file, load_config, resolve_path
from npstructfeat.pore import load_residue_features
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。

    参数说明
    ----
    --config:
        YAML 配置文件路径。

    --thresholds:
        要测试的径向距离阈值列表，单位 Å。
        默认测试 [15, 18, 20, 22, 25]。
    """

    parser = argparse.ArgumentParser(
        description="Analyze sensitivity of inner residue selection to radius threshold."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="config/default.yaml",
        help="Path to YAML config file. Default: config/default.yaml",
    )

    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=[15.0, 18.0, 20.0, 22.0, 25.0],
        help="Radius thresholds to test, in Angstrom. Default: 15 18 20 22 25",
    )

    return parser.parse_args()


def build_features_for_threshold(
    residue_geo_df: pd.DataFrame,
    config: Dict[str, Any],
    threshold: float,
) -> Dict[str, Any]:
    """
    对单个阈值进行候选孔道内壁残基筛选，并聚合结构特征。

    参数
    ----
    residue_geo_df : pandas.DataFrame
        已经添加 radial_distance、theta、z_norm 等几何特征的残基表。

    config : dict
        YAML 配置字典。

    threshold : float
        当前测试的径向距离阈值，单位 Å。

    返回
    ----
    features : dict
        当前阈值下的 nanopore-level 结构统计结果。

    说明
    ----
    该函数不会保存 inner_candidate_residues.csv。
    它只在内存中筛选并统计，避免覆盖主流程输出结果。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")

    axis_mode = config.get("pore", {}).get("axis_mode", "z_axis")
    center_mode = config.get("pore", {}).get("center_mode", "xy_mean")

    if axis_mode != "z_axis":
        raise ValueError(
            f"当前敏感性分析仅支持 axis_mode='z_axis'，但得到: {axis_mode}"
        )

    # 根据当前阈值筛选候选孔道内壁残基
    inner_df = residue_geo_df[
        residue_geo_df["radial_distance"] <= threshold
    ].copy()

    total_residue_count = len(residue_geo_df)
    inner_count = len(inner_df)

    # 如果某个阈值筛不到任何残基，返回空值，避免后续统计报错
    if inner_count == 0:
        return {
            "nanopore_id": nanopore_id,
            "pdb_id": pdb_id,
            "axis_mode": axis_mode,
            "center_mode": center_mode,
            "threshold_A": threshold,
            "total_residue_count": total_residue_count,
            "inner_candidate_residue_count": 0,
            "inner_candidate_residue_ratio": 0.0,
            "chain_count": 0,
            "inner_candidate_count_per_chain_mean": None,
            "inner_candidate_count_per_chain_std": None,
            "inner_candidate_count_per_chain_min": None,
            "inner_candidate_count_per_chain_max": None,
            "z_min_A": None,
            "z_max_A": None,
            "pore_region_length_approx_A": None,
            "radial_distance_min_A": None,
            "radial_distance_max_A": None,
            "radial_distance_mean_A": None,
            "radial_distance_std_A": None,
            "inner_net_charge": None,
            "inner_mean_charge": None,
            "inner_positive_residue_count": None,
            "inner_negative_residue_count": None,
            "inner_neutral_residue_count": None,
            "inner_mean_hydrophobicity": None,
            "inner_std_hydrophobicity": None,
            "inner_aromatic_count": None,
            "inner_aromatic_ratio": None,
            "inner_polar_count": None,
            "inner_polar_ratio": None,
            "inner_residue_type_count": None,
        }

    # 每条链候选残基数量，用于检查七聚体对称性
    chain_counts = inner_df.groupby("chain_id")["residue_number"].count()

    z_min = float(inner_df["z"].min())
    z_max = float(inner_df["z"].max())
    z_range = z_max - z_min

    radial_min = float(inner_df["radial_distance"].min())
    radial_max = float(inner_df["radial_distance"].max())
    radial_mean = float(inner_df["radial_distance"].mean())
    radial_std = float(inner_df["radial_distance"].std())

    inner_net_charge = float(inner_df["charge"].sum())
    inner_mean_charge = float(inner_df["charge"].mean())

    inner_mean_hydrophobicity = float(inner_df["hydrophobicity"].mean())
    inner_std_hydrophobicity = float(inner_df["hydrophobicity"].std())

    aromatic_count = int(inner_df["is_aromatic"].sum())
    aromatic_ratio = float(inner_df["is_aromatic"].mean())

    polar_count = int(inner_df["is_polar"].sum())
    polar_ratio = float(inner_df["is_polar"].mean())

    positive_count = int((inner_df["charge"] > 0).sum())
    negative_count = int((inner_df["charge"] < 0).sum())
    neutral_count = int((inner_df["charge"] == 0).sum())

    residue_type_count = int(inner_df["residue_name"].nunique())

    # center_x / center_y 在 add_basic_geometry_features 中已经写入每一行
    center_x = float(inner_df["center_x"].iloc[0])
    center_y = float(inner_df["center_y"].iloc[0])

    features = {
        "nanopore_id": nanopore_id,
        "pdb_id": pdb_id,
        "axis_mode": axis_mode,
        "center_mode": center_mode,
        "center_x": center_x,
        "center_y": center_y,

        # 当前测试阈值
        "threshold_A": threshold,

        # 数量特征
        "total_residue_count": total_residue_count,
        "inner_candidate_residue_count": inner_count,
        "inner_candidate_residue_ratio": inner_count / total_residue_count,
        "chain_count": int(inner_df["chain_id"].nunique()),

        # 各链分布特征
        "inner_candidate_count_per_chain_mean": float(chain_counts.mean()),
        "inner_candidate_count_per_chain_std": float(chain_counts.std()),
        "inner_candidate_count_per_chain_min": int(chain_counts.min()),
        "inner_candidate_count_per_chain_max": int(chain_counts.max()),

        # z 方向覆盖范围
        "z_min_A": z_min,
        "z_max_A": z_max,
        "pore_region_length_approx_A": z_range,

        # 径向距离统计
        "radial_distance_min_A": radial_min,
        "radial_distance_max_A": radial_max,
        "radial_distance_mean_A": radial_mean,
        "radial_distance_std_A": radial_std,

        # 电荷特征
        "inner_net_charge": inner_net_charge,
        "inner_mean_charge": inner_mean_charge,
        "inner_positive_residue_count": positive_count,
        "inner_negative_residue_count": negative_count,
        "inner_neutral_residue_count": neutral_count,

        # 疏水性特征
        "inner_mean_hydrophobicity": inner_mean_hydrophobicity,
        "inner_std_hydrophobicity": inner_std_hydrophobicity,

        # 芳香族/极性特征
        "inner_aromatic_count": aromatic_count,
        "inner_aromatic_ratio": aromatic_ratio,
        "inner_polar_count": polar_count,
        "inner_polar_ratio": polar_ratio,

        # 残基组成复杂度
        "inner_residue_type_count": residue_type_count,
    }

    return features


def run_threshold_sensitivity(
    config: Dict[str, Any],
    thresholds: List[float],
) -> pd.DataFrame:
    """
    对多个阈值运行敏感性分析。

    参数
    ----
    config : dict
        YAML 配置字典。

    thresholds : list[float]
        需要测试的阈值列表，单位 Å。

    返回
    ----
    sensitivity_df : pandas.DataFrame
        阈值敏感性分析结果表。
        每一行对应一个 threshold。
    """

    pore_config = config.get("pore", {})
    center_mode = pore_config.get("center_mode", "xy_mean")

    # 读取 residue_features.csv
    residue_df = load_residue_features(config)

    # 只需要对原始残基表做一次几何特征添加
    # 后续不同阈值都复用这个 residue_geo_df
    residue_geo_df = add_basic_geometry_features(
        residue_df=residue_df,
        center_mode=center_mode,
    )

    rows = []

    for threshold in thresholds:
        features = build_features_for_threshold(
            residue_geo_df=residue_geo_df,
            config=config,
            threshold=float(threshold),
        )
        rows.append(features)

    sensitivity_df = pd.DataFrame(rows)

    # 按阈值从小到大排序
    sensitivity_df = sensitivity_df.sort_values("threshold_A").reset_index(drop=True)

    return sensitivity_df


def get_threshold_sensitivity_output_file(config: Dict[str, Any]) -> Path:
    """
    生成阈值敏感性分析结果的输出路径。

    默认输出到：
        data/processed/nanopore_features/<PDB_ID>_threshold_sensitivity.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_file : Path
        输出 CSV 文件路径。
    """

    pdb_id = config["input"]["pdb_id"].upper()

    output_dir = resolve_path(config["output"]["nanopore_feature_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{pdb_id}_threshold_sensitivity.csv"

    return output_file


def print_sensitivity_summary(sensitivity_df: pd.DataFrame) -> None:
    """
    打印敏感性分析的核心结果。

    参数
    ----
    sensitivity_df : pandas.DataFrame
        阈值敏感性分析表。
    """

    display_cols = [
        "threshold_A",
        "inner_candidate_residue_count",
        "inner_candidate_residue_ratio",
        "inner_net_charge",
        "inner_mean_hydrophobicity",
        "inner_aromatic_ratio",
        "inner_polar_ratio",
        "pore_region_length_approx_A",
        "inner_candidate_count_per_chain_min",
        "inner_candidate_count_per_chain_max",
    ]

    print("\n阈值敏感性分析核心结果：")

    # 使用 to_string()，避免 pandas 在终端中省略中间列
    print(
        sensitivity_df[display_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}"
        )
    )

    print("\n解释建议：")
    print(
        "1. inner_candidate_residue_count 应随 threshold 增大而增加；\n"
        "2. inner_net_charge、inner_mean_hydrophobicity 等特征若变化过剧烈，"
        "说明阈值选择对结构特征影响较大；\n"
        "3. 每条链候选残基数量的 min/max 如果差距很小，说明七聚体对称性较好；\n"
        "4. 当前分析仍属于几何粗筛验证，不等价于 HOLE/CHAP 的严格孔道识别。"
    )


def main():
    """
    主函数。

    执行流程：
    1. 读取配置；
    2. 检查必要字段；
    3. 打印配置摘要；
    4. 对多个 threshold 运行敏感性分析；
    5. 打印核心结果；
    6. 保存 CSV。
    """

    args = parse_args()

    config = load_config(args.config)

    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    print_config_summary(config)

    thresholds = args.thresholds

    print("\n待测试阈值 threshold_A:")
    print(thresholds)

    sensitivity_df = run_threshold_sensitivity(
        config=config,
        thresholds=thresholds,
    )

    print_sensitivity_summary(sensitivity_df)

    output_file = get_threshold_sensitivity_output_file(config)
    sensitivity_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(f"\n阈值敏感性分析结果已保存到：{output_file}")


if __name__ == "__main__":
    main()