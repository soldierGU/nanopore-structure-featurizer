"""
03_select_inner_residues.py

本脚本用于筛选候选孔道内壁残基。

运行方式：
    python scripts/03_select_inner_residues.py

或者指定配置文件：
    python scripts/03_select_inner_residues.py --config config/default.yaml

输入文件：
    data/processed/residue_features/5JZT_residue_features.csv

输出文件：
    data/processed/inner_residues/5JZT_inner_candidate_residues.csv

当前筛选方法：
    radial_distance <= inner_radius_threshold

其中 radial_distance 是残基 Cα 到孔道中心轴的距离。
第一版假设孔道中心轴沿 z 轴。
"""

import argparse

from npstructfeat.io import load_config
from npstructfeat.pore import (
    save_inner_candidate_residues,
    select_inner_candidate_residues,
    summarize_inner_candidates,
)
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Select candidate pore-lining residues."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="config/default.yaml",
        help="Path to YAML config file. Default: config/default.yaml",
    )

    return parser.parse_args()


def main():
    """
    主函数。

    执行流程：
    1. 读取 YAML 配置；
    2. 检查配置字段；
    3. 打印配置摘要；
    4. 筛选候选孔道内壁残基；
    5. 打印统计结果；
    6. 保存 CSV 文件。
    """

    args = parse_args()

    # 1. 读取配置
    config = load_config(args.config)

    # 2. 检查必要字段
    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    # 3. 打印配置摘要
    print_config_summary(config)

    # 4. 筛选候选孔道内壁残基
    inner_df = select_inner_candidate_residues(config)

    print("\n候选孔道内壁残基表预览：")
    print(inner_df.head())

    print("\n候选孔道内壁残基总体统计：")
    print(f"候选残基总数: {len(inner_df)}")
    print(f"涉及链数量: {inner_df['chain_id'].nunique()}")
    print(f"残基类型数量: {inner_df['residue_name'].nunique()}")

    print("\n径向距离统计 radial_distance:")
    print(inner_df["radial_distance"].describe())

    print("\nz 坐标范围:")
    print(f"z_min = {inner_df['z'].min():.3f}")
    print(f"z_max = {inner_df['z'].max():.3f}")
    print(f"z_range = {inner_df['z'].max() - inner_df['z'].min():.3f} Å")

    print("\n各链候选孔道内壁残基统计：")
    summary_df = summarize_inner_candidates(inner_df)
    print(summary_df)

    # 5. 保存结果
    output_file = save_inner_candidate_residues(config)

    print(f"\n候选孔道内壁残基表已保存到：{output_file}")

    print(
        "\n提醒：当前结果是基于 radial_distance 阈值的几何粗筛，"
        "不是 HOLE/CHAP 等专业工具识别的严格 pore-lining residues。"
    )


if __name__ == "__main__":
    main()