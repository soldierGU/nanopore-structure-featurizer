"""
io.py

本文件负责工程中的输入输出操作，主要包括：

1. 读取 YAML 配置文件；
2. 检查输入结构文件是否存在；
3. 创建输出目录；
4. 根据配置生成标准化输出路径。

为什么单独写 io.py？
因为后面的每个脚本都会用到这些功能，例如：
- 解析结构时需要读取 config/default.yaml；
- 保存 residue_features.csv 前需要确认输出目录存在；
- 运行 pipeline 前需要检查 raw CIF 文件是否存在。

把这些公共逻辑放在 io.py，可以避免在多个脚本里重复写代码。
"""

from pathlib import Path
from typing import Any, Dict

import yaml


def load_config(config_path: str | Path) -> Dict[str, Any]:
    """
    读取 YAML 配置文件。

    参数
    ----
    config_path : str | Path
        YAML 配置文件路径，例如：
        "config/default.yaml"

    返回
    ----
    config : dict
        读取后的配置字典。

    示例
    ----
    config = load_config("config/default.yaml")
    pdb_id = config["input"]["pdb_id"]

    注意
    ----
    yaml.safe_load() 比 yaml.load() 更安全，
    一般读取普通配置文件时优先使用 safe_load。
    """

    config_path = Path(config_path)

    # 如果是相对路径，转换为相对于项目根目录的绝对路径
    if not config_path.is_absolute():
        project_root = get_project_root()
        config_path = project_root / config_path

    if not config_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_path}")

    if not config_path.is_file():
        raise ValueError(f"配置路径不是文件: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ValueError(f"配置文件为空: {config_path}")

    return config


def get_project_root() -> Path:
    """
    获取工程根目录。

    当前假设：
    这个文件位于：
        src/npstructfeat/io.py

    那么工程根目录是：
        src/npstructfeat/io.py
        -> npstructfeat
        -> src
        -> project_root

    返回
    ----
    project_root : Path
        工程根目录的绝对路径。

    为什么需要这个函数？
    因为我们希望无论从哪里运行脚本，
    都能正确找到 data/、config/、outputs/ 等目录。
    """

    current_file = Path(__file__).resolve()

    # parents[0] = src/npstructfeat
    # parents[1] = src
    # parents[2] = 工程根目录
    project_root = current_file.parents[2]

    return project_root


def resolve_path(path: str | Path, project_root: str | Path | None = None) -> Path:
    """
    将配置文件中的相对路径转换为绝对路径。

    参数
    ----
    path : str | Path
        文件或目录路径。
        例如：
        "data/raw/5JZT/5JZT_assembly.cif"

    project_root : str | Path | None
        工程根目录。
        如果为 None，则自动调用 get_project_root()。

    返回
    ----
    resolved_path : Path
        绝对路径。

    说明
    ----
    YAML 文件中通常写相对路径，这样工程更容易迁移。
    但 Python 实际读写文件时，使用绝对路径更稳妥。
    """

    path = Path(path)

    # 如果本来就是绝对路径，直接返回
    if path.is_absolute():
        return path

    if project_root is None:
        project_root = get_project_root()

    project_root = Path(project_root)

    return project_root / path


def check_input_structure_file(config: Dict[str, Any]) -> Path:
    """
    检查输入的 PDB/mmCIF 结构文件是否存在。

    参数
    ----
    config : dict
        从 YAML 文件读取到的配置字典。

    返回
    ----
    structure_file : Path
        结构文件的绝对路径。

    会检查的内容
    ----
    1. config 中是否存在 input 字段；
    2. input 中是否存在 structure_file 字段；
    3. structure_file 对应的文件是否真实存在；
    4. 该路径是否确实是文件，而不是目录。
    """

    if "input" not in config:
        raise KeyError("配置文件中缺少 input 字段。")

    if "structure_file" not in config["input"]:
        raise KeyError("配置文件 input 中缺少 structure_file 字段。")

    structure_file = resolve_path(config["input"]["structure_file"])

    if not structure_file.exists():
        raise FileNotFoundError(
            f"输入结构文件不存在: {structure_file}\n"
            "请检查 config/default.yaml 中的 input.structure_file 是否正确。"
        )

    if not structure_file.is_file():
        raise ValueError(f"输入结构路径不是文件: {structure_file}")

    return structure_file


def ensure_dir(dir_path: str | Path) -> Path:
    """
    确保某个目录存在。

    如果目录不存在，则自动创建。
    如果目录已经存在，则什么也不做。

    参数
    ----
    dir_path : str | Path
        需要创建或检查的目录路径。

    返回
    ----
    dir_path : Path
        目录的 Path 对象。
    """

    dir_path = Path(dir_path)

    # parents=True 表示可以递归创建多级目录
    # exist_ok=True 表示目录已存在时不报错
    dir_path.mkdir(parents=True, exist_ok=True)

    return dir_path


def prepare_output_dirs(config: Dict[str, Any]) -> Dict[str, Path]:
    """
    根据 YAML 配置创建输出目录。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_dirs : dict[str, Path]
        输出目录字典，方便后续代码直接使用。

    示例
    ----
    output_dirs = prepare_output_dirs(config)
    residue_dir = output_dirs["residue_feature_dir"]
    """

    if "output" not in config:
        raise KeyError("配置文件中缺少 output 字段。")

    output_config = config["output"]

    output_dirs: Dict[str, Path] = {}

    for key, path_str in output_config.items():
        # YAML 中写的是相对路径，这里转换成绝对路径
        abs_path = resolve_path(path_str)

        # 自动创建目录
        ensure_dir(abs_path)

        output_dirs[key] = abs_path

    return output_dirs


def get_output_file(
    config: Dict[str, Any],
    output_dir_key: str,
    suffix: str,
) -> Path:
    """
    根据 PDB ID 和输出目录类型，生成标准化输出文件路径。

    参数
    ----
    config : dict
        YAML 配置字典。

    output_dir_key : str
        输出目录字段名。
        例如：
        "residue_feature_dir"
        "chain_summary_dir"

    suffix : str
        文件后缀名部分。
        例如：
        "residue_features.csv"
        "chain_summary.csv"

    返回
    ----
    output_file : Path
        标准化输出路径。

    示例
    ----
    get_output_file(
        config,
        output_dir_key="residue_feature_dir",
        suffix="residue_features.csv"
    )

    如果 pdb_id = 5JZT，则输出：
    data/processed/residue_features/5JZT_residue_features.csv
    """

    pdb_id = config["input"]["pdb_id"].upper()

    if "output" not in config:
        raise KeyError("配置文件中缺少 output 字段。")

    if output_dir_key not in config["output"]:
        raise KeyError(f"配置文件 output 中缺少 {output_dir_key} 字段。")

    output_dir = resolve_path(config["output"][output_dir_key])
    ensure_dir(output_dir)

    output_file = output_dir / f"{pdb_id}_{suffix}"

    return output_file