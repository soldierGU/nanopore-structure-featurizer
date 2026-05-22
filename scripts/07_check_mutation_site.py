"""
07_check_mutation_site.py

本脚本用于检查突变位点是否适合后续 FoldX 突变建模。

当前默认用于检查：
    T232K

主要检查内容：
1. 5JZT 中 A-G 七条链的 residue_number=232 是否存在；
2. 该位点实际残基是否为 THR；
3. 该位点是否位于 candidate inner residues 中；
4. 该位点的 radial_distance、z、z_norm 等空间信息。

运行方式：
    python scripts/07_check_mutation_site.py

指定突变：
    python scripts/07_check_mutation_site.py --mutation T232K

指定配置文件：
    python scripts/07_check_mutation_site.py --config config/default.yaml --mutation T232K

输出：
    data/processed/mutation_sites/5JZT_T232K_site_check.csv
"""

import argparse

from npstructfeat.io import load_config
from npstructfeat.mutation import (
    check_mutation_set_sites,
    save_mutation_set_site_check,
    summarize_mutation_set_site_check,
)
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。

    参数
    ----
    --config:
        YAML 配置文件路径。

    --mutation:
        突变字符串，例如 T232K、K238Q。
        当前只支持简单单点突变格式。
    """

    parser = argparse.ArgumentParser(
        description="Check mutation site before FoldX mutant modeling."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="config/default.yaml",
        help="Path to YAML config file. Default: config/default.yaml",
    )

    parser.add_argument(
        "--mutation",
        type=str,
        default="T232K",
        help="Mutation ID, e.g., T232K. Default: T232K",
    )

    return parser.parse_args()


def main():
    """
    主函数。

    执行流程：
    1. 读取 YAML 配置；
    2. 检查必要字段；
    3. 打印配置摘要；
    4. 检查突变位点；
    5. 打印检查表；
    6. 打印简要统计；
    7. 保存 CSV 文件。
    """

    args = parse_args()

    config = load_config(args.config)

    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    print_config_summary(config)

    mutation_id = args.mutation.strip().upper()

    print(f"\n开始检查突变位点: {mutation_id}")

    site_df = check_mutation_set_sites(
        config=config,
        mutation_id=mutation_id,
    )

    print("\n突变位点检查结果：")
    display_cols = [
        "mutation_id",
        "chain_id",
        "residue_number",
        "expected_residue",
        "observed_residue",
        "target_residue",
        "is_site_found",
        "is_expected_match",
        "is_inner_candidate",
        "radial_distance",
        "z",
        "z_norm",
    ]

    print(
        site_df[display_cols].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    summary = summarize_mutation_set_site_check(site_df)

    print("\n检查统计：")
    for key, value in summary.items():
        print(f"{key}: {value}")

    output_file = save_mutation_set_site_check(
        config=config,
        mutation_id=mutation_id,
    )

    print(f"\n突变位点检查表已保存到：{output_file}")

    print("\n判断建议：")

    if summary["all_sites_found"] and summary["all_expected_match"]:
        print(
            f"1. {mutation_id} 的 WT 位点在所有检查链中均存在，"
            "且实际残基与预期残基一致。"
        )
    else:
        print(
            "1. 位点检查未完全通过。请检查 residue_number 是否与 PDB 编号一致，"
            "必要时需要做序列编号映射。"
        )

    if summary["all_inner_candidate"]:
        print(
            "2. 该突变位点在所有链中均属于 candidate inner residues，"
            "说明它位于当前几何筛选的候选孔道区域内。"
        )
    elif summary["inner_candidate_count"] > 0:
        print(
            "2. 该突变位点只在部分链中属于 candidate inner residues。"
            "需要检查结构对称性或筛选阈值。"
        )
    else:
        print(
            "2. 该突变位点未出现在 candidate inner residues 中。"
            "这不一定说明突变无意义，但需要重新评估阈值或位点位置。"
        )

    print(
        "\n注意：通过该检查并不代表已经完成突变建模。"
        "它只是 FoldX BuildModel 之前的编号和位置验证步骤。"
    )


if __name__ == "__main__":
    main()
