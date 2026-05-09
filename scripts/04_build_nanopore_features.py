"""
04_build_nanopore_features.py

本脚本用于将候选孔道内壁残基表聚合为 nanopore-level 结构特征表。

运行方式：
    python scripts/04_build_nanopore_features.py

或者指定配置文件：
    python scripts/04_build_nanopore_features.py --config config/default.yaml

输入文件：
    data/processed/inner_residues/5JZT_inner_candidate_residues.csv

输出文件：
    data/processed/nanopore_features/5JZT_nanopore_structure_features.csv

该输出文件后续可以与实验样本表拼接，用于：
    待测物 + 实验条件 + 纳米孔结构特征 → 传感响应回归
"""

import argparse

from npstructfeat.features import (
    build_nanopore_structure_features,
    save_nanopore_structure_features,
)
from npstructfeat.io import load_config
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Build nanopore-level structural feature table."
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
    4. 聚合 nanopore-level 结构特征；
    5. 打印结果；
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

    # 4. 构建纳米孔整体结构特征
    feature_df = build_nanopore_structure_features(config)

    print("\n纳米孔整体结构特征表：")
    print(feature_df.T)

    # 5. 保存结果
    output_file = save_nanopore_structure_features(config)

    print(f"\n纳米孔整体结构特征表已保存到：{output_file}")

    print(
        "\n说明：该表基于 candidate inner residues 聚合得到，"
        "当前适合作为结构先验特征表。"
    )


if __name__ == "__main__":
    main()