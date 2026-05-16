"""
10_compare_wt_mutant_features.py

本脚本用于比较 WT 与 FoldX mutant 的 nanopore-level 结构特征差异。

当前默认比较：
    WT:    5JZT_nanopore_structure_features.csv
    Mutant: 5JZT_T232K_nanopore_structure_features.csv

输出：
    data/processed/delta_features/5JZT_WT_vs_T232K_delta_features.csv

运行方式：
    python scripts/10_compare_wt_mutant_features.py

也可以指定输入文件：
    python scripts/10_compare_wt_mutant_features.py \
        --wt-feature data/processed/nanopore_features/5JZT_nanopore_structure_features.csv \
        --mutant-feature data/processed/nanopore_features/5JZT_T232K_nanopore_structure_features.csv
"""

import argparse
from pathlib import Path

from npstructfeat.compare import (
    compare_wt_mutant_feature_files,
    save_delta_features,
)
from npstructfeat.io import resolve_path


DEFAULT_FEATURES_TO_COMPARE = [
    # 残基数量和候选区域
    "total_residue_count",
    "inner_candidate_residue_count",
    "inner_candidate_residue_ratio",
    "chain_count",
    "inner_candidate_count_per_chain_mean",
    "inner_candidate_count_per_chain_std",
    "inner_candidate_count_per_chain_min",
    "inner_candidate_count_per_chain_max",

    # z 方向覆盖范围
    "z_min_A",
    "z_max_A",
    "pore_region_length_approx_A",

    # 径向距离统计
    "radial_distance_min_A",
    "radial_distance_max_A",
    "radial_distance_mean_A",
    "radial_distance_std_A",

    # 电荷相关
    "inner_net_charge",
    "inner_mean_charge",
    "inner_positive_residue_count",
    "inner_negative_residue_count",
    "inner_neutral_residue_count",

    # 疏水性相关
    "inner_mean_hydrophobicity",
    "inner_std_hydrophobicity",

    # 芳香族/极性
    "inner_aromatic_count",
    "inner_aromatic_ratio",
    "inner_polar_count",
    "inner_polar_ratio",

    # 残基类型复杂度
    "inner_residue_type_count",
]


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Compare WT and mutant nanopore structural features."
    )

    parser.add_argument(
        "--wt-feature",
        type=str,
        default="data/processed/nanopore_features/5JZT_nanopore_structure_features.csv",
        help="WT nanopore structure feature table.",
    )

    parser.add_argument(
        "--mutant-feature",
        type=str,
        default=(
            "data/processed/nanopore_features/"
            "5JZT_T232K_nanopore_structure_features.csv"
        ),
        help="Mutant nanopore structure feature table.",
    )

    parser.add_argument(
        "--wt-id",
        type=str,
        default="aerolysin_WT",
        help="WT nanopore ID.",
    )

    parser.add_argument(
        "--mutant-id",
        type=str,
        default="aerolysin_T232K",
        help="Mutant nanopore ID.",
    )

    parser.add_argument(
        "--template-pdb-id",
        type=str,
        default="5JZT",
        help="Template PDB ID.",
    )

    parser.add_argument(
        "--mutation",
        type=str,
        default="T232K",
        help="Mutation ID.",
    )

    parser.add_argument(
        "--modeling-method",
        type=str,
        default="FoldX",
        help="Mutant modeling method.",
    )

    parser.add_argument(
        "--compare-all-numeric",
        action="store_true",
        help=(
            "Compare all shared numeric features instead of the predefined "
            "feature list."
        ),
    )

    return parser.parse_args()


def print_key_delta_summary(delta_df):
    """
    打印关键 delta features。
    """

    key_features = [
        "inner_net_charge",
        "inner_positive_residue_count",
        "inner_negative_residue_count",
        "inner_neutral_residue_count",
        "inner_mean_charge",
        "inner_mean_hydrophobicity",
        "inner_polar_ratio",
        "inner_aromatic_ratio",
        "inner_candidate_residue_count",
        "radial_distance_mean_A",
        "pore_region_length_approx_A",
    ]

    key_df = delta_df[delta_df["feature_name"].isin(key_features)].copy()

    # 保持 key_features 的显示顺序
    key_df["feature_order"] = key_df["feature_name"].apply(
        lambda x: key_features.index(x) if x in key_features else 999
    )
    key_df = key_df.sort_values("feature_order").drop(columns=["feature_order"])

    display_cols = [
        "feature_name",
        "wt_value",
        "mutant_value",
        "delta_value",
        "relative_delta",
    ]

    print("\n关键 WT-mutant 差异特征：")
    print(
        key_df[display_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.6f}",
        )
    )


def main():
    """
    主函数。

    执行流程：
    1. 读取 WT feature table；
    2. 读取 mutant feature table；
    3. 计算 mutant - WT；
    4. 打印关键差异；
    5. 保存 delta feature table。
    """

    args = parse_args()

    wt_feature_file = resolve_path(args.wt_feature)
    mutant_feature_file = resolve_path(args.mutant_feature)

    print("WT feature file:")
    print(wt_feature_file)

    print("\nMutant feature file:")
    print(mutant_feature_file)

    if not wt_feature_file.exists():
        raise FileNotFoundError(f"WT feature file 不存在: {wt_feature_file}")

    if not mutant_feature_file.exists():
        raise FileNotFoundError(f"Mutant feature file 不存在: {mutant_feature_file}")

    if args.compare_all_numeric:
        feature_names = None
        print("\n比较模式：所有共同数值型特征")
    else:
        feature_names = DEFAULT_FEATURES_TO_COMPARE
        print("\n比较模式：预定义核心结构特征")

    delta_df = compare_wt_mutant_feature_files(
        wt_feature_file=wt_feature_file,
        mutant_feature_file=mutant_feature_file,
        wt_id=args.wt_id,
        mutant_id=args.mutant_id,
        template_pdb_id=args.template_pdb_id,
        mutation_id=args.mutation,
        modeling_method=args.modeling_method,
        feature_names=feature_names,
    )

    print_key_delta_summary(delta_df)

    output_file = save_delta_features(
        delta_df=delta_df,
        template_pdb_id=args.template_pdb_id,
        mutation_id=args.mutation,
    )

    print(f"\nWT-mutant delta features 已保存到：{output_file}")

    print("\n解释提醒：")
    print(
        "1. delta_value = mutant_value - wt_value。\n"
        "2. 对 T232K，最应关注 inner_net_charge、"
        "inner_positive_residue_count、inner_mean_hydrophobicity 和 inner_polar_ratio。\n"
        "3. 如果 inner_candidate_residue_count 也发生变化，说明突变后几何筛选集合发生了变化。"
    )


if __name__ == "__main__":
    main()