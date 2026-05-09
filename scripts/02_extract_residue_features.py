"""
02_extract_residue_features.py

本脚本用于从 PDB/mmCIF 结构文件中提取残基级特征表。

运行方式：
    python scripts/02_extract_residue_features.py

或者指定配置文件：
    python scripts/02_extract_residue_features.py --config config/default.yaml

输出文件：
    data/processed/residue_features/5JZT_residue_features.csv

该文件后续用于：
1. 筛选候选孔道内壁残基；
2. 构建残基接触图；
3. 聚合纳米孔整体结构特征；
4. 与待测物和实验条件特征融合。
"""

import argparse

from npstructfeat.features import extract_residue_features, save_residue_features
from npstructfeat.io import load_config
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Extract residue-level features from nanopore structure."
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
    2. 检查必要配置字段；
    3. 打印配置摘要；
    4. 提取残基级特征；
    5. 打印结果预览；
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

    # 4. 提取残基特征
    residue_df = extract_residue_features(config)

    print("\n残基级特征表预览：")
    print(residue_df.head())

    print("\n残基级特征表统计：")
    print(f"总残基数: {len(residue_df)}")
    print(f"链数量: {residue_df['chain_id'].nunique()}")
    print(f"残基类型数量: {residue_df['residue_name'].nunique()}")
    print(f"缺失 CA 残基数: {(residue_df['has_ca'] == 0).sum()}")

    print("\n各链残基数量：")
    print(residue_df.groupby("chain_id")["residue_number"].count())

    # 5. 保存结果
    output_file = save_residue_features(config)

    print(f"\n残基级特征表已保存到：{output_file}")


if __name__ == "__main__":
    main()