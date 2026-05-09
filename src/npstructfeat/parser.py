"""
parser.py

本文件负责解析 PDB/mmCIF 结构文件，并提取基础结构统计信息。

当前阶段的主要功能：
1. 根据配置文件读取 mmCIF/PDB 结构；
2. 统计每条链的残基数量；
3. 统计每条链的原子数量；
4. 统计标准氨基酸残基数量；
5. 统计缺失 CA 原子的残基数量；
6. 输出 chain-level summary 表格。

为什么先做 chain summary？
因为在真正提取纳米孔结构特征之前，必须先确认：
- 下载的是不是 biological assembly；
- 是否包含完整七聚体；
- 每条链是否大致完整；
- 结构文件是否能被 Biopython 正常解析。

如果这一步有问题，后续孔道内壁筛选和结构特征提取都会不可靠。
"""

from pathlib import Path
from typing import Any, Dict, Iterable

import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from npstructfeat.io import check_input_structure_file, get_output_file


# 20 种标准氨基酸三字母缩写
# 注意：PDB/mmCIF 里残基名称通常是三字母大写，例如 ALA、GLY、LYS。
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}


def load_structure(config: Dict[str, Any]) -> Structure:
    """
    根据配置文件读取结构文件。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    structure : Bio.PDB.Structure.Structure
        Biopython 解析后的结构对象。

    支持格式
    ----
    1. cif / mmcif
    2. pdb

    说明
    ----
    你的 5JZT Biological Assembly 解压后是 .cif 文件，
    所以当前主要会使用 MMCIFParser。
    """

    structure_file = check_input_structure_file(config)

    pdb_id = config["input"]["pdb_id"].upper()
    file_format = config["input"].get("file_format", "cif").lower()

    if file_format in {"cif", "mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif file_format == "pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            f"不支持的结构文件格式: {file_format}\n"
            "当前仅支持 cif/mmcif/pdb。"
        )

    structure = parser.get_structure(pdb_id, structure_file)

    return structure


def get_selected_model(structure: Structure, model_index: int = 0):
    """
    从结构中选择一个 model。

    参数
    ----
    structure : Bio.PDB.Structure.Structure
        Biopython 结构对象。

    model_index : int
        选择第几个 model。
        通常 X-ray 和 cryo-EM 结构只有一个 model，即 index=0。
        NMR 结构可能有多个 model。

    返回
    ----
    model : Bio.PDB.Model.Model
        被选中的模型。

    为什么要有这个函数？
    因为 PDB/mmCIF 文件的层级结构通常是：
        Structure -> Model -> Chain -> Residue -> Atom

    我们现在先只处理第一个 model。
    """

    models = list(structure.get_models())

    if len(models) == 0:
        raise ValueError("结构文件中没有任何 model。")

    if model_index < 0 or model_index >= len(models):
        raise IndexError(
            f"model_index={model_index} 超出范围。"
            f"该结构共有 {len(models)} 个 model。"
        )

    return models[model_index]


def is_standard_residue(residue: Residue) -> bool:
    """
    判断一个 Residue 是否为标准氨基酸残基。

    参数
    ----
    residue : Bio.PDB.Residue.Residue
        Biopython 的残基对象。

    返回
    ----
    bool
        True 表示标准氨基酸残基；
        False 表示水分子、配体、离子、非标准残基等。

    说明
    ----
    residue.id[0] == " " 通常表示普通蛋白残基。
    但为了更稳妥，这里还检查 residue.resname 是否属于 20 种标准氨基酸。
    """

    # residue.id[0] 不是空格时，通常表示 hetero residue，例如水、配体、离子
    if residue.id[0] != " ":
        return False

    return residue.resname.strip().upper() in STANDARD_AMINO_ACIDS


def iter_residues(chain: Chain) -> Iterable[Residue]:
    """
    遍历一条链中的所有残基。

    参数
    ----
    chain : Bio.PDB.Chain.Chain
        Biopython 的链对象。

    返回
    ----
    Iterable[Residue]
        残基迭代器。

    这个函数目前比较简单，主要是为了让代码语义更清晰。
    """

    for residue in chain:
        yield residue


def summarize_chain(chain: Chain) -> Dict[str, Any]:
    """
    统计单条链的信息。

    参数
    ----
    chain : Bio.PDB.Chain.Chain
        Biopython 链对象。

    返回
    ----
    summary : dict
        当前链的统计信息。

    统计内容
    ----
    chain_id:
        链 ID，例如 A、B、C。

    total_residue_count:
        链中所有残基数量，包括标准残基、非标准残基、水、配体等。

    standard_residue_count:
        标准氨基酸残基数量。

    nonstandard_residue_count:
        非标准残基数量。

    atom_count:
        当前链中的总原子数量。

    ca_missing_count:
        标准氨基酸残基中缺失 CA 原子的数量。

    first_residue_number / last_residue_number:
        标准氨基酸残基的编号范围。
    """

    total_residue_count = 0
    standard_residue_count = 0
    nonstandard_residue_count = 0
    atom_count = 0
    ca_missing_count = 0

    standard_residue_numbers = []

    for residue in iter_residues(chain):
        total_residue_count += 1
        atom_count += len(list(residue.get_atoms()))

        if is_standard_residue(residue):
            standard_residue_count += 1
            standard_residue_numbers.append(residue.id[1])

            if "CA" not in residue:
                ca_missing_count += 1
        else:
            nonstandard_residue_count += 1

    if standard_residue_numbers:
        first_residue_number = min(standard_residue_numbers)
        last_residue_number = max(standard_residue_numbers)
    else:
        first_residue_number = None
        last_residue_number = None

    summary = {
        "chain_id": chain.id,
        "total_residue_count": total_residue_count,
        "standard_residue_count": standard_residue_count,
        "nonstandard_residue_count": nonstandard_residue_count,
        "atom_count": atom_count,
        "ca_missing_count": ca_missing_count,
        "first_residue_number": first_residue_number,
        "last_residue_number": last_residue_number,
    }

    return summary


def build_chain_summary(config: Dict[str, Any]) -> pd.DataFrame:
    """
    构建链级别结构统计表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    df : pandas.DataFrame
        链级别统计表。

    表格每一行代表一条链。
    """

    structure = load_structure(config)

    model_index = config.get("structure", {}).get("use_model_index", 0)
    model = get_selected_model(structure, model_index=model_index)

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")

    rows = []

    for chain in model:
        summary = summarize_chain(chain)

        # 增加 PDB 和纳米孔信息，方便以后合并多个结构
        summary["pdb_id"] = pdb_id
        summary["nanopore_id"] = nanopore_id
        summary["model_index"] = model_index

        rows.append(summary)

    df = pd.DataFrame(rows)

    # 调整列顺序，让表格更清晰
    columns = [
        "pdb_id",
        "nanopore_id",
        "model_index",
        "chain_id",
        "total_residue_count",
        "standard_residue_count",
        "nonstandard_residue_count",
        "atom_count",
        "ca_missing_count",
        "first_residue_number",
        "last_residue_number",
    ]

    df = df[columns]

    return df


def save_chain_summary(config: Dict[str, Any]) -> Path:
    """
    构建并保存链级别结构统计表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    output_file : Path
        保存后的 CSV 文件路径。
    """

    df = build_chain_summary(config)

    output_file = get_output_file(
        config=config,
        output_dir_key="chain_summary_dir",
        suffix="chain_summary.csv",
    )

    df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file