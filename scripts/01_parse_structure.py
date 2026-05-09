"""
01_parse_structure.py

本脚本用于解析 PDB/mmCIF 结构文件，并生成链级别统计表。

运行方式：
    python scripts/01_parse_structure.py

或者指定配置文件：
    python scripts/01_parse_structure.py --config config/default.yaml

输出文件：
    data/processed/chain_summary/5JZT_chain_summary.csv

该文件用于检查：
1. Biological Assembly 是否包含完整七聚体；
2. 每条链的残基数量是否一致；
3. 是否存在明显缺失 CA 原子的残基；
4. 是否存在大量非标准残基。
"""

import argparse

from npstructfeat.io import load_config
from npstructfeat.parser import build_chain_summary, save_chain_summary
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。

    当前只支持一个参数：
    --config

    如果不指定，则默认使用：
    config/default.yaml
    """

    parser = argparse.ArgumentParser(
        description="Parse nanopore structure and build chain summary."
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
    1. 读取配置文件；
    2. 检查必要配置字段；
    3. 打印配置摘要；
    4. 解析结构；
    5. 构建链级别 summary；
    6. 保存 CSV 文件；
    7. 在终端打印结果预览。
    """

    args = parse_args()

    # 1. 读取 YAML 配置
    config = load_config(args.config)

    # 2. 检查配置字段
    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    # 3. 打印配置摘要，便于确认当前处理的是哪个 PDB 文件
    print_config_summary(config)

    # 4. 构建链级别统计表
    chain_summary_df = build_chain_summary(config)

    print("\n链级别统计结果：")
    print(chain_summary_df)

    # 5. 保存结果
    output_file = save_chain_summary(config)

    print(f"\n链级别统计表已保存到：{output_file}")

    # 6. 根据链数量给出一个简单提示
    chain_count = len(chain_summary_df)
    print(f"\n检测到链数量：{chain_count}")

    if chain_count == 7:
        print("链数量为 7，符合 aerolysin 七聚体孔道的预期。")
    else:
        print(
            "链数量不是 7。请检查下载的是否为 Biological Assembly，"
            "或者该结构本身是否不是七聚体。"
        )


if __name__ == "__main__":
    main()