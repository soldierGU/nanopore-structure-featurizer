"""
00_check_config.py

这个脚本用于检查：

1. YAML 配置文件是否能成功读取；
2. 输入的 CIF/PDB 文件是否存在；
3. 输出目录是否能自动创建；
4. 当前配置是否符合预期。

运行方式：
    python scripts/00_check_config.py

或者指定配置文件：
    python scripts/00_check_config.py --config config/default.yaml
"""
import src
import argparse

from src.npstructfeat.io import (
    check_input_structure_file,
    load_config,
    prepare_output_dirs,
)
from src.npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。

    argparse 是 Python 标准库，用于从命令行读取参数。

    这里我们只定义一个参数：
    --config

    如果用户没有指定，就默认使用：
    config/default.yaml
    """

    parser = argparse.ArgumentParser(
        description="Check YAML config and project paths."
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

    脚本运行逻辑：
    1. 读取命令行参数；
    2. 读取 YAML 配置；
    3. 检查必要字段；
    4. 打印配置摘要；
    5. 检查输入结构文件；
    6. 创建输出目录。
    """

    args = parse_args()

    # 1. 读取 YAML 配置文件
    config = load_config(args.config)

    # 2. 检查配置文件中是否包含必要的一级字段
    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    # 3. 打印配置摘要，方便人工检查
    print_config_summary(config)

    # 4. 检查输入 CIF/PDB 文件是否存在
    structure_file = check_input_structure_file(config)
    print(f"\n输入结构文件检查通过: {structure_file}")

    # 5. 创建输出目录
    output_dirs = prepare_output_dirs(config)
    print("\n输出目录已准备完成：")
    for key, path in output_dirs.items():
        print(f"{key}: {path}")

    print("\n配置检查完成。")


if __name__ == "__main__":
    main()