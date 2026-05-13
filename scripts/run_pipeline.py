"""
run_pipeline.py

一键运行 Nanopore Structure Featurizer 的完整流程。

运行方式：
    python scripts/run_pipeline.py

或者指定配置文件：
    python scripts/run_pipeline.py --config config/default.yaml

执行步骤：
1. 检查配置和输入结构文件；
2. 解析结构并生成链级别统计表；
3. 提取残基级特征表；
4. 筛选候选孔道内壁残基；
5. 聚合纳米孔整体结构特征表。

这个脚本适合以后处理新的 PDB 文件，例如：
    python scripts/run_pipeline.py --config config/9FM6.yaml
"""

import argparse

from npstructfeat.io import (
    check_input_structure_file,
    load_config,
    prepare_output_dirs,
)
from npstructfeat.parser import save_chain_summary
from npstructfeat.features import (
    save_residue_features,
    save_nanopore_structure_features,
)
from npstructfeat.pore import save_inner_candidate_residues
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Run full nanopore structure featurization pipeline."
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

    一次性完成从结构文件到结构特征表的完整流程。
    """

    args = parse_args()

    # 1. 读取配置
    config = load_config(args.config)

    # 2. 检查必要配置字段
    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    # 3. 打印配置摘要
    print_config_summary(config)

    # 4. 检查输入结构文件
    structure_file = check_input_structure_file(config)
    print(f"\n[OK] 输入结构文件存在: {structure_file}")

    # 5. 创建输出目录
    prepare_output_dirs(config)
    print("[OK] 输出目录已准备完成")

    # 6. 解析结构，保存 chain summary
    print("\n[Step 1/4] 解析结构并生成链级别统计表...")
    chain_summary_file = save_chain_summary(config)
    print(f"[OK] 链级别统计表已保存: {chain_summary_file}")

    # 7. 提取 residue-level features
    print("\n[Step 2/4] 提取残基级结构特征...")
    residue_feature_file = save_residue_features(config)
    print(f"[OK] 残基级特征表已保存: {residue_feature_file}")

    # 8. 筛选 candidate inner residues
    print("\n[Step 3/4] 筛选候选孔道内壁残基...")
    inner_residue_file = save_inner_candidate_residues(config)
    print(f"[OK] 候选孔道内壁残基表已保存: {inner_residue_file}")

    # 9. 聚合 nanopore-level features
    print("\n[Step 4/4] 构建纳米孔整体结构特征表...")
    nanopore_feature_file = save_nanopore_structure_features(config)
    print(f"[OK] 纳米孔结构特征表已保存: {nanopore_feature_file}")

    print("\n完整流程运行结束。")


if __name__ == "__main__":
    main()

# 使用示例 当前工程下terminal或者powershell运行run_pipeline.py --config config/9FM6.yaml
# 其中9FM6为纳米孔Pdb编号，运行时请将9FM6替换为对应编号