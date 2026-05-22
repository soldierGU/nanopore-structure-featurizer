"""
foldx.py

本模块用于准备 FoldX 相关输入文件和辅助检查。

当前阶段主要功能：
1. 将 mmCIF/PDB 结构转换为 FoldX 输入用 PDB 文件；
2. 只保留标准氨基酸残基；
3. 可选择保留指定链；
4. 输出 FoldX 输入准备报告；
5. 检查突变位点在导出的 PDB 中是否仍然存在且编号正确。

注意：
FoldX 是外部命令行程序，不是 Python 包。
本模块当前只负责准备输入文件，不直接运行 FoldX RepairPDB 或 BuildModel。
"""

import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from npstructfeat.io import check_input_structure_file, resolve_path
from npstructfeat.parser import get_selected_model
from npstructfeat.residue_props import STANDARD_AMINO_ACIDS
from npstructfeat.mutation import normalize_mutation_label, parse_mutation_set

class FoldXProteinSelect(Select):
    """
    Biopython PDBIO 的结构选择器。

    作用：
    在保存 PDB 文件时，只保留我们希望 FoldX 处理的部分。

    当前规则：
    1. 只保留指定链；
    2. 只保留标准氨基酸残基；
    3. 去除水分子、离子、小分子、配体等非标准残基。

    为什么要这样做？
    FoldX 对输入 PDB 比较敏感。
    对当前点突变建模任务来说，我们只需要蛋白七聚体结构。
    """

    def __init__(self, chains_to_keep: Optional[List[str]] = None):
        """
        参数
        ----
        chains_to_keep : list[str] | None
            需要保留的链 ID。
            如果为 None，则保留所有链。
        """

        super().__init__()

        if chains_to_keep is None:
            self.chains_to_keep = None
        else:
            self.chains_to_keep = set(str(chain_id) for chain_id in chains_to_keep)

    def accept_chain(self, chain: Chain) -> int:
        """
        判断是否保留某条链。

        返回 1 表示保留；
        返回 0 表示丢弃。
        """

        if self.chains_to_keep is None:
            return 1

        return int(chain.id in self.chains_to_keep)

    def accept_residue(self, residue: Residue) -> int:
        """
        判断是否保留某个残基。

        只保留标准氨基酸残基。

        Biopython 中：
        residue.id[0] == " " 通常表示普通蛋白残基；
        其他情况通常表示水、配体、离子、非标准残基等。
        """

        if residue.id[0] != " ":
            return 0

        resname = residue.resname.strip().upper()

        return int(resname in STANDARD_AMINO_ACIDS)


def load_structure_for_foldx(config: Dict[str, Any]) -> Structure:
    """
    根据配置文件读取结构文件。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    structure : Bio.PDB.Structure.Structure
        Biopython 结构对象。
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


def get_foldx_input_dir(config: Dict[str, Any]) -> Path:
    """
    获取 FoldX 输入目录。

    默认使用：
        data/foldx/input/<PDB_ID>_WT/

    如果 config 中存在 foldx.input_dir，则使用该路径作为根目录。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        当前结构对应的 FoldX input 目录。
    """

    pdb_id = config["input"]["pdb_id"].upper()

    foldx_config = config.get("foldx", {})
    input_root = foldx_config.get("input_dir", "data/foldx/input")

    input_dir = resolve_path(input_root) / f"{pdb_id}_WT"
    input_dir.mkdir(parents=True, exist_ok=True)

    return input_dir


def export_foldx_input_pdb(
    config: Dict[str, Any],
    chains_to_keep: Optional[List[str]] = None,
) -> Path:
    """
    将当前结构导出为 FoldX 输入用 PDB 文件。

    参数
    ----
    config : dict
        YAML 配置字典。

    chains_to_keep : list[str] | None
        需要保留的链。
        如果为 None，则保留结构中的所有链。

    返回
    ----
    output_pdb : Path
        导出的 PDB 文件路径。

    输出示例
    ----
    data/foldx/input/5JZT_WT/5JZT_assembly.pdb
    """

    structure = load_structure_for_foldx(config)

    model_index = config.get("structure", {}).get("use_model_index", 0)
    model = get_selected_model(structure, model_index=model_index)

    pdb_id = config["input"]["pdb_id"].upper()
    output_dir = get_foldx_input_dir(config)

    output_pdb = output_dir / f"{pdb_id}_assembly.pdb"

    io = PDBIO()
    io.set_structure(model)

    selector = FoldXProteinSelect(chains_to_keep=chains_to_keep)
    io.save(str(output_pdb), select=selector)

    return output_pdb


def summarize_exported_pdb(
    pdb_file: str | Path,
    pdb_id: str,
    nanopore_id: str,
) -> pd.DataFrame:
    """
    读取导出的 PDB 文件，并统计链级别信息。

    参数
    ----
    pdb_file : str | Path
        导出的 FoldX 输入 PDB 文件。

    pdb_id : str
        PDB ID，例如 5JZT。

    nanopore_id : str
        纳米孔 ID，例如 aerolysin_WT。

    返回
    ----
    summary_df : pandas.DataFrame
        每条链一行的统计表。

    统计内容
    ----
    - chain_id
    - standard_residue_count
    - atom_count
    - first_residue_number
    - last_residue_number
    """

    pdb_file = Path(pdb_file)

    if not pdb_file.exists():
        raise FileNotFoundError(f"导出的 PDB 文件不存在: {pdb_file}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)

    model = get_selected_model(structure, model_index=0)

    rows = []

    for chain in model:
        standard_residue_numbers = []
        atom_count = 0

        for residue in chain:
            atom_count += len(list(residue.get_atoms()))

            if residue.id[0] == " " and residue.resname.strip().upper() in STANDARD_AMINO_ACIDS:
                standard_residue_numbers.append(residue.id[1])

        if standard_residue_numbers:
            first_residue_number = min(standard_residue_numbers)
            last_residue_number = max(standard_residue_numbers)
        else:
            first_residue_number = None
            last_residue_number = None

        rows.append(
            {
                "pdb_id": pdb_id,
                "nanopore_id": nanopore_id,
                "chain_id": chain.id,
                "standard_residue_count": len(standard_residue_numbers),
                "atom_count": atom_count,
                "first_residue_number": first_residue_number,
                "last_residue_number": last_residue_number,
            }
        )

    summary_df = pd.DataFrame(rows)

    return summary_df


def check_residue_in_exported_pdb(
    pdb_file: str | Path,
    residue_number: int,
    expected_residue: str,
) -> pd.DataFrame:
    """
    检查导出的 PDB 文件中指定编号残基是否存在。

    参数
    ----
    pdb_file : str | Path
        导出的 FoldX 输入 PDB 文件。

    residue_number : int
        需要检查的残基编号，例如 232。

    expected_residue : str
        预期三字母残基名，例如 THR。

    返回
    ----
    site_df : pandas.DataFrame
        每条链一行，记录该位点是否存在以及是否匹配预期残基。
    """

    pdb_file = Path(pdb_file)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("exported_pdb", pdb_file)

    model = get_selected_model(structure, model_index=0)

    rows = []

    expected_residue = expected_residue.strip().upper()

    for chain in model:
        hit = None

        for residue in chain:
            if residue.id[0] != " ":
                continue

            if residue.id[1] == residue_number:
                hit = residue
                break

        if hit is None:
            rows.append(
                {
                    "chain_id": chain.id,
                    "residue_number": residue_number,
                    "expected_residue": expected_residue,
                    "observed_residue": None,
                    "is_site_found": False,
                    "is_expected_match": False,
                }
            )
        else:
            observed_residue = hit.resname.strip().upper()

            rows.append(
                {
                    "chain_id": chain.id,
                    "residue_number": residue_number,
                    "expected_residue": expected_residue,
                    "observed_residue": observed_residue,
                    "is_site_found": True,
                    "is_expected_match": observed_residue == expected_residue,
                }
            )

    site_df = pd.DataFrame(rows)

    return site_df


def prepare_foldx_input_report(
    config: Dict[str, Any],
    output_pdb: str | Path,
    mutation_id: str = "T232K",
) -> pd.DataFrame:
    """
    ?? FoldX ???????

    ??????????????
    ???????????????????????????
    """

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")
    mutation_infos = parse_mutation_set(mutation_id)

    summary_df = summarize_exported_pdb(
        pdb_file=output_pdb,
        pdb_id=pdb_id,
        nanopore_id=nanopore_id,
    )

    report_frames = []

    for mutation_info in mutation_infos:
        site_df = check_residue_in_exported_pdb(
            pdb_file=output_pdb,
            residue_number=mutation_info["residue_number"],
            expected_residue=mutation_info["original_aa_3"],
        )

        # ??????????????????????????
        report_df = summary_df.merge(
            site_df,
            on="chain_id",
            how="left",
        )
        report_df["mutation_set_id"] = normalize_mutation_label(mutation_id)
        report_df["mutation_id"] = mutation_info["mutation_id"]
        report_df["target_residue"] = mutation_info["target_aa_3"]
        report_df["output_pdb"] = str(output_pdb)
        report_frames.append(report_df)

    report_df = pd.concat(report_frames, ignore_index=True)

    columns = [
        "pdb_id",
        "nanopore_id",
        "mutation_set_id",
        "mutation_id",
        "chain_id",
        "standard_residue_count",
        "atom_count",
        "first_residue_number",
        "last_residue_number",
        "residue_number",
        "expected_residue",
        "observed_residue",
        "target_residue",
        "is_site_found",
        "is_expected_match",
        "output_pdb",
    ]

    return report_df[columns]


def save_foldx_input_report(
    config: Dict[str, Any],
    output_pdb: str | Path,
    mutation_id: str = "T232K",
) -> Path:
    """
    保存 FoldX 输入准备报告。

    输出路径：
        data/foldx/input/<PDB_ID>_WT/<PDB_ID>_foldx_input_report.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    output_pdb : str | Path
        导出的 FoldX 输入 PDB 文件。

    mutation_id : str
        突变 ID。

    返回
    ----
    Path
        保存后的报告路径。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    input_dir = get_foldx_input_dir(config)

    report_df = prepare_foldx_input_report(
        config=config,
        output_pdb=output_pdb,
        mutation_id=mutation_id,
    )

    report_file = input_dir / f"{pdb_id}_foldx_input_report.csv"
    report_df.to_csv(report_file, index=False, encoding="utf-8-sig")

    return report_file

def get_foldx_executable(config: Dict[str, Any]) -> Path:
    """
    获取 FoldX 可执行文件路径。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        FoldX 可执行程序路径。

    要求
    ----
    config/default.yaml 中建议包含：

    foldx:
      executable: tools/foldx/foldx_1_20270131.exe
    """

    foldx_config = config.get("foldx", {})
    executable = foldx_config.get("executable", None)

    if executable is None:
        raise KeyError(
            "配置文件中缺少 foldx.executable 字段。\n"
            "请在 config/default.yaml 中添加：\n"
            "foldx:\n"
            "  executable: tools/foldx/foldx_1_20270131.exe"
        )

    executable_path = resolve_path(executable)

    if not executable_path.exists():
        raise FileNotFoundError(f"FoldX 可执行文件不存在: {executable_path}")

    return executable_path


def get_foldx_repaired_pdb(config: Dict[str, Any]) -> Path:
    """
    获取 RepairPDB 生成的 repaired PDB 文件路径。

    默认预期路径：
        data/foldx/input/<PDB_ID>_WT/<PDB_ID>_assembly_Repair.pdb

    例如：
        data/foldx/input/5JZT_WT/5JZT_assembly_Repair.pdb

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        repaired PDB 文件路径。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    input_dir = get_foldx_input_dir(config)

    repaired_pdb = input_dir / f"{pdb_id}_assembly_Repair.pdb"

    if not repaired_pdb.exists():
        raise FileNotFoundError(
            f"FoldX repaired PDB 不存在: {repaired_pdb}\n"
            "请先完成 FoldX RepairPDB：\n"
            f"foldx --command=RepairPDB --pdb={pdb_id}_assembly.pdb"
        )

    return repaired_pdb


def get_foldx_mutant_dir(
    config: Dict[str, Any],
    mutation_id: str,
) -> Path:
    """
    获取 FoldX mutant 工作目录。

    默认路径：
        data/foldx/mutants/<PDB_ID>_<MUTATION_ID>/

    例如：
        data/foldx/mutants/5JZT_T232K/

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    返回
    ----
    Path
        FoldX mutant 工作目录。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    mutation_label = normalize_mutation_label(mutation_id)

    foldx_config = config.get("foldx", {})
    mutant_root = foldx_config.get("mutant_dir", "data/foldx/mutants")

    # 用规范化后的标签命名目录，避免 "/" 之类字符污染路径。
    mutant_dir = resolve_path(mutant_root) / f"{pdb_id}_{mutation_label}"
    mutant_dir.mkdir(parents=True, exist_ok=True)

    return mutant_dir


def get_modeled_mutant_dir(
    config: Dict[str, Any],
    mutation_id: str,
) -> Path:
    """
    获取最终整理后的 modeled mutant 输出目录。

    默认路径：
        data/modeled/<PDB_ID>_<MUTATION_ID>/

    例如：
        data/modeled/5JZT_T232K/

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    返回
    ----
    Path
        modeled mutant 输出目录。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    mutation_label = normalize_mutation_label(mutation_id)

    foldx_config = config.get("foldx", {})
    modeled_root = foldx_config.get("modeled_dir", "data/modeled")

    modeled_dir = resolve_path(modeled_root) / f"{pdb_id}_{mutation_label}"
    modeled_dir.mkdir(parents=True, exist_ok=True)

    return modeled_dir


def write_foldx_individual_list(
    config: Dict[str, Any],
    mutation_id: str,
    chains: List[str],
) -> Path:
    """
    生成 FoldX BuildModel 所需的 individual_list.txt 文件。

    支持单点和多点突变。

    示例 1：
        T232K + A-G
        -> TA232K,TB232K,TC232K,TD232K,TE232K,TF232K,TG232K;

    示例 2：
        T232K/K238Q + A-G
        -> TA232K,TB232K,...,TG232K,KA238Q,KB238Q,...,KG238Q;
    """

    mutation_infos = parse_mutation_set(mutation_id)

    mutant_dir = get_foldx_mutant_dir(
        config=config,
        mutation_id=mutation_id,
    )

    mutation_items = []

    for mutation_info in mutation_infos:
        original_aa = mutation_info["original_aa_1"]
        target_aa = mutation_info["target_aa_1"]
        residue_number = mutation_info["residue_number"]

        for chain_id in chains:
            chain_id = str(chain_id)
            mutation_items.append(
                f"{original_aa}{chain_id}{residue_number}{target_aa}"
            )

    mutation_line = ",".join(mutation_items) + ";"

    individual_list_file = mutant_dir / "individual_list.txt"
    individual_list_file.write_text(mutation_line + "\n", encoding="utf-8")

    return individual_list_file


def prepare_foldx_buildmodel_workspace(
    config: Dict[str, Any],
    mutation_id: str,
    chains: List[str],
    overwrite_individual_list: bool = False,
) -> Dict[str, Path]:
    """
    准备 FoldX BuildModel 工作目录。

    主要操作：
    1. 检查 repaired PDB 是否存在；
    2. 将 repaired PDB 复制到 mutant 工作目录；
    3. 检查或生成 individual_list.txt。

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    chains : list[str]
        需要突变的链。

    overwrite_individual_list : bool
        如果 True，则覆盖已有 individual_list.txt；
        如果 False，已有文件会被保留。

    返回
    ----
    dict
        包含工作目录和关键文件路径。
    """

    repaired_pdb = get_foldx_repaired_pdb(config)

    mutant_dir = get_foldx_mutant_dir(
        config=config,
        mutation_id=mutation_id,
    )

    # 将 repaired PDB 复制到 mutant 工作目录
    repaired_pdb_in_mutant_dir = mutant_dir / repaired_pdb.name
    shutil.copy2(repaired_pdb, repaired_pdb_in_mutant_dir)

    individual_list_file = mutant_dir / "individual_list.txt"

    if overwrite_individual_list or not individual_list_file.exists():
        individual_list_file = write_foldx_individual_list(
            config=config,
            mutation_id=mutation_id,
            chains=chains,
        )

    return {
        "mutant_dir": mutant_dir,
        "repaired_pdb": repaired_pdb,
        "repaired_pdb_in_mutant_dir": repaired_pdb_in_mutant_dir,
        "individual_list_file": individual_list_file,
    }


def run_foldx_buildmodel(
    config: Dict[str, Any],
    mutation_id: str,
    chains: List[str],
    overwrite_individual_list: bool = False,
) -> Dict[str, Any]:
    """
    运行 FoldX BuildModel。

    注意：
    该函数会真正调用 FoldX，因此只有在 RepairPDB 完成后才能运行。

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    chains : list[str]
        需要突变的链。

    overwrite_individual_list : bool
        是否覆盖已有 individual_list.txt。

    返回
    ----
    result : dict
        包含命令、返回码、stdout、stderr、工作目录等信息。
    """

    foldx_executable = get_foldx_executable(config)

    workspace = prepare_foldx_buildmodel_workspace(
        config=config,
        mutation_id=mutation_id,
        chains=chains,
        overwrite_individual_list=overwrite_individual_list,
    )

    mutant_dir = workspace["mutant_dir"]
    repaired_pdb_in_mutant_dir = workspace["repaired_pdb_in_mutant_dir"]
    individual_list_file = workspace["individual_list_file"]

    cmd = [
        str(foldx_executable),
        "--command=BuildModel",
        f"--pdb={repaired_pdb_in_mutant_dir.name}",
        f"--mutant-file={individual_list_file.name}",
    ]

    completed = subprocess.run(
        cmd,
        cwd=str(mutant_dir),
        capture_output=True,
        text=True,
    )

    return {
        "cmd": cmd,
        "cwd": mutant_dir,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "workspace": workspace,
    }


def find_foldx_mutant_pdb(
    config: Dict[str, Any],
    mutation_id: str,
) -> Path:
    """
    查找 FoldX BuildModel 生成的 mutant PDB 文件。

    FoldX 通常会输出类似：
        5JZT_assembly_Repair_1.pdb

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    返回
    ----
    Path
        找到的 mutant PDB 文件。

    注意
    ----
    如果有多个候选文件，会优先选择包含 '_1.pdb' 且不是 repaired 输入文件的 PDB。
    """

    mutant_dir = get_foldx_mutant_dir(
        config=config,
        mutation_id=mutation_id,
    )

    pdb_files = sorted(mutant_dir.glob("*.pdb"))

    if len(pdb_files) == 0:
        raise FileNotFoundError(f"FoldX mutant 目录中没有找到 PDB 文件: {mutant_dir}")

    # 排除 copied repaired PDB
    candidate_files = [
        p for p in pdb_files
        if not p.name.endswith("_Repair.pdb")
    ]

    # 优先选择 FoldX 常见输出 _1.pdb
    preferred = [
        p for p in candidate_files
        if p.stem.endswith("_1")
    ]

    if len(preferred) == 1:
        return preferred[0]

    if len(preferred) > 1:
        return preferred[0]

    if len(candidate_files) == 1:
        return candidate_files[0]

    raise FileNotFoundError(
        f"无法唯一确定 FoldX mutant PDB。\n"
        f"目录: {mutant_dir}\n"
        f"候选文件: {[p.name for p in candidate_files]}"
    )


def collect_foldx_mutant_model(
    config: Dict[str, Any],
    mutation_id: str,
) -> Path:
    """
    将 FoldX 生成的 mutant PDB 整理到 data/modeled 目录。

    输出路径：
        data/modeled/<PDB_ID>_<MUTATION_ID>/<PDB_ID>_<MUTATION_ID>_model.pdb

    例如：
        data/modeled/5JZT_T232K/5JZT_T232K_model.pdb

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变 ID，例如 T232K。

    返回
    ----
    Path
        整理后的 mutant model PDB 路径。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    mutation_label = normalize_mutation_label(mutation_id)

    foldx_mutant_pdb = find_foldx_mutant_pdb(
        config=config,
        mutation_id=mutation_id,
    )

    modeled_dir = get_modeled_mutant_dir(
        config=config,
        mutation_id=mutation_id,
    )

    # 输出文件名同样必须使用规范化标签，否则多突变字符串会被当成子路径。
    output_model = modeled_dir / f"{pdb_id}_{mutation_label}_model.pdb"

    shutil.copy2(foldx_mutant_pdb, output_model)

    return output_model


def validate_mutation_set_in_pdb(
    pdb_file: str | Path,
    mutation_id: str,
    chains: List[str],
) -> pd.DataFrame:
    """
    检查 mutant PDB 中多个突变位点是否全部变成目标残基。

    对 T232K/K238Q，会检查：
        A-G:232 是否为 LYS
        A-G:238 是否为 GLN
    """

    mutation_infos = parse_mutation_set(mutation_id)

    pdb_file = Path(pdb_file)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mutant_model", pdb_file)
    model = get_selected_model(structure, model_index=0)

    chain_dict = {chain.id: chain for chain in model}

    rows = []

    for mutation_info in mutation_infos:
        single_mutation_id = mutation_info["mutation_id"]
        residue_number = mutation_info["residue_number"]
        target_residue = mutation_info["target_aa_3"]

        for chain_id in chains:
            chain_id = str(chain_id)

            if chain_id not in chain_dict:
                rows.append(
                    {
                        "mutation_id": single_mutation_id,
                        "chain_id": chain_id,
                        "residue_number": residue_number,
                        "target_residue": target_residue,
                        "observed_residue": None,
                        "is_chain_found": False,
                        "is_site_found": False,
                        "is_target_match": False,
                    }
                )
                continue

            chain = chain_dict[chain_id]
            hit = None

            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if residue.id[1] == residue_number:
                    hit = residue
                    break

            if hit is None:
                rows.append(
                    {
                        "mutation_id": single_mutation_id,
                        "chain_id": chain_id,
                        "residue_number": residue_number,
                        "target_residue": target_residue,
                        "observed_residue": None,
                        "is_chain_found": True,
                        "is_site_found": False,
                        "is_target_match": False,
                    }
                )
            else:
                observed_residue = hit.resname.strip().upper()

                rows.append(
                    {
                        "mutation_id": single_mutation_id,
                        "chain_id": chain_id,
                        "residue_number": residue_number,
                        "target_residue": target_residue,
                        "observed_residue": observed_residue,
                        "is_chain_found": True,
                        "is_site_found": True,
                        "is_target_match": observed_residue == target_residue,
                    }
                )

    return pd.DataFrame(rows)
