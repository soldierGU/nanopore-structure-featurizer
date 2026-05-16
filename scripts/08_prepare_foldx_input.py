"""
08_prepare_foldx_input.py

本脚本用于为 FoldX 准备输入 PDB 文件。

当前任务：
1. 读取 config/default.yaml 中的 5JZT_assembly.cif；
2. 将 mmCIF 转换成 PDB；
3. 只保留标准氨基酸残基；
4. 保留 A-G 七条链；
5. 输出 FoldX 输入用 PDB 文件；
6. 生成 FoldX input preparation report；
7. 检查 T232 位点在导出 PDB 中是否仍然为 THR。

运行方式：
    python scripts/08_prepare_foldx_input.py

指定突变：
    python scripts/08_prepare_foldx_input.py --mutation T232K

输出：
    data/foldx/input/5JZT_WT/5JZT_assembly.pdb
    data/foldx/input/5JZT_WT/5JZT_foldx_input_report.csv
"""

import argparse

from npstructfeat.foldx import (
    export_foldx_input_pdb,
    prepare_foldx_input_report,
    save_foldx_input_report,
)
from npstructfeat.io import load_config
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Prepare FoldX input PDB file from nanopore structure."
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
        help="Mutation ID for site check, e.g., T232K. Default: T232K",
    )

    parser.add_argument(
        "--chains",
        type=str,
        nargs="+",
        default=["A", "B", "C", "D", "E", "F", "G"],
        help="Chains to keep in FoldX input PDB. Default: A B C D E F G",
    )

    return parser.parse_args()


def main():
    """
    主函数。

    执行流程：
    1. 读取配置；
    2. 检查配置字段；
    3. 打印配置摘要；
    4. 导出 FoldX 输入 PDB；
    5. 生成并打印 FoldX 输入准备报告；
    6. 保存报告。
    """

    args = parse_args()

    config = load_config(args.config)

    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    print_config_summary(config)

    mutation_id = args.mutation.strip().upper()
    chains_to_keep = args.chains

    print("\n开始准备 FoldX 输入文件")
    print(f"Mutation site check: {mutation_id}")
    print(f"Chains to keep: {chains_to_keep}")

    output_pdb = export_foldx_input_pdb(
        config=config,
        chains_to_keep=chains_to_keep,
    )

    print(f"\nFoldX 输入 PDB 已导出：{output_pdb}")

    report_df = prepare_foldx_input_report(
        config=config,
        output_pdb=output_pdb,
        mutation_id=mutation_id,
    )

    print("\nFoldX 输入准备报告：")
    display_cols = [
        "chain_id",
        "standard_residue_count",
        "first_residue_number",
        "last_residue_number",
        "residue_number",
        "expected_residue",
        "observed_residue",
        "target_residue",
        "is_site_found",
        "is_expected_match",
    ]

    print(
        report_df[display_cols].to_string(
            index=False,
        )
    )

    report_file = save_foldx_input_report(
        config=config,
        output_pdb=output_pdb,
        mutation_id=mutation_id,
    )

    print(f"\nFoldX 输入准备报告已保存到：{report_file}")

    all_site_found = bool(report_df["is_site_found"].all())
    all_expected_match = bool(report_df["is_expected_match"].all())

    print("\n检查结论：")

    if all_site_found and all_expected_match:
        print(
            f"1. 导出的 PDB 中，所有保留链的 {mutation_id} WT 位点均存在，"
            "且残基与预期一致。"
        )
        print("2. 该 PDB 可以作为 FoldX RepairPDB 的输入。")
    else:
        print(
            "1. 导出的 PDB 中存在位点缺失或残基不匹配问题。\n"
            "2. 请先检查链 ID、残基编号或 CIF→PDB 转换结果，暂不建议进入 FoldX RepairPDB。"
        )

    print(
        "\n下一步：\n"
        "在 PowerShell 中手动运行 FoldX RepairPDB，确认 FoldX 可以正常处理该 PDB。"
    )


if __name__ == "__main__":
    main()