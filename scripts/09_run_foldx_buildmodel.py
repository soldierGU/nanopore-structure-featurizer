"""
09_run_foldx_buildmodel.py

本脚本用于运行 FoldX BuildModel，构建突变体结构模型。

当前默认任务：
    5JZT WT repaired structure
    → T232K 七聚体突变模型

前提条件：
1. 已运行 scripts/08_prepare_foldx_input.py；
2. 已手动运行 FoldX RepairPDB；
3. 已生成：
   data/foldx/input/5JZT_WT/5JZT_assembly_Repair.pdb

运行方式：
    python scripts/09_run_foldx_buildmodel.py

指定突变：
    python scripts/09_run_foldx_buildmodel.py --mutation T232K

输出：
    data/foldx/mutants/5JZT_T232K/
    data/modeled/5JZT_T232K/5JZT_T232K_model.pdb

注意：
本脚本会真正调用 FoldX BuildModel。
只有在 RepairPDB 完成后才能运行。
"""

import argparse
from pathlib import Path

from npstructfeat.foldx import (
    collect_foldx_mutant_model,
    get_foldx_mutant_dir,
    prepare_foldx_buildmodel_workspace,
    run_foldx_buildmodel,
    validate_mutation_set_in_pdb,
)
from npstructfeat.io import load_config
from npstructfeat.utils import print_config_summary, require_config_keys
from npstructfeat.mutation import normalize_mutation_label

def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Run FoldX BuildModel for nanopore mutant modeling."
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

    parser.add_argument(
        "--chains",
        type=str,
        nargs="+",
        default=["A", "B", "C", "D", "E", "F", "G"],
        help="Chains to mutate. Default: A B C D E F G",
    )

    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help=(
            "Only prepare BuildModel workspace and individual_list.txt, "
            "do not run FoldX."
        ),
    )

    parser.add_argument(
        "--overwrite-individual-list",
        action="store_true",
        help="Overwrite existing individual_list.txt.",
    )

    return parser.parse_args()


def main():
    """
    主函数。

    执行流程：
    1. 读取配置；
    2. 检查 repaired PDB 是否存在；
    3. 准备 FoldX BuildModel 工作目录；
    4. 如果 --prepare-only，则只准备文件，不运行 FoldX；
    5. 否则调用 FoldX BuildModel；
    6. 收集 FoldX 输出 mutant PDB；
    7. 检查 A-G 七条链是否成功突变为目标残基。
    """

    args = parse_args()

    config = load_config(args.config)

    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    print_config_summary(config)

    mutation_id = args.mutation.strip().upper()
    chains = args.chains

    print("\nFoldX BuildModel 设置：")
    print(f"Mutation: {mutation_id}")
    print(f"Chains: {chains}")

    print("\n准备 FoldX BuildModel 工作目录...")

    workspace = prepare_foldx_buildmodel_workspace(
        config=config,
        mutation_id=mutation_id,
        chains=chains,
        overwrite_individual_list=args.overwrite_individual_list,
    )

    print(f"Mutant work dir: {workspace['mutant_dir']}")
    print(f"Repaired PDB copied to: {workspace['repaired_pdb_in_mutant_dir']}")
    print(f"Individual list: {workspace['individual_list_file']}")

    individual_text = Path(workspace["individual_list_file"]).read_text(
        encoding="utf-8"
    ).strip()

    print(f"\nindividual_list.txt 内容：\n{individual_text}")

    if args.prepare_only:
        print("\n当前为 --prepare-only 模式，不运行 FoldX BuildModel。")
        print("RepairPDB 完成后，可去掉 --prepare-only 重新运行本脚本。")
        return

    print("\n开始运行 FoldX BuildModel...")

    result = run_foldx_buildmodel(
        config=config,
        mutation_id=mutation_id,
        chains=chains,
        overwrite_individual_list=False,
    )

    print("\nFoldX command:")
    print(" ".join(result["cmd"]))

    print("\nFoldX return code:")
    print(result["returncode"])

    print("\nFoldX stdout:")
    print(result["stdout"])

    if result["stderr"]:
        print("\nFoldX stderr:")
        print(result["stderr"])

    if result["returncode"] != 0:
        print("\nFoldX BuildModel 返回非 0 状态。请检查 stdout/stderr。")
        return

    print("\n收集 FoldX mutant PDB...")

    mutant_model = collect_foldx_mutant_model(
        config=config,
        mutation_id=mutation_id,
    )

    print(f"整理后的 mutant model 已保存到：{mutant_model}")

    print("\n验证 mutant 位点是否突变成功...")

    validation_df = validate_mutation_set_in_pdb(
        pdb_file=mutant_model,
        mutation_id=mutation_id,
        chains=chains,
    )

    print(
        validation_df.to_string(
            index=False,
        )
    )

    mutant_dir = get_foldx_mutant_dir(
        config=config,
        mutation_id=mutation_id,
    )

    mutation_label = normalize_mutation_label(mutation_id)
    validation_file = mutant_dir / f"{mutation_label}_mutant_validation.csv"
    validation_df.to_csv(validation_file, index=False, encoding="utf-8-sig")

    print(f"\n突变验证表已保存到：{validation_file}")

    if bool(validation_df["is_target_match"].all()):
        print("\n所有指定链的目标位点均已突变为目标残基。")
        print("下一步可以为 mutant model 建立独立 YAML，并运行现有 pipeline。")
    else:
        print("\n存在未成功突变的链或位点。请先检查 FoldX 输出和 individual_list.txt。")


if __name__ == "__main__":
    main()