"""
mutation.py

本模块用于处理纳米孔蛋白突变位点相关分析。

当前主要功能：
1. 检查指定突变位点在 WT 结构中是否存在；
2. 检查该位点的实际残基是否与预期残基一致；
3. 检查该位点是否属于 candidate inner residues；
4. 输出 mutation site check 表格，为后续 FoldX 突变建模做准备。

为什么需要这个模块？
在进行 T232K、K238Q 等突变体建模之前，必须确认：

- PDB 文件中的 residue_number 是否与文献编号一致；
- 七聚体中每条链对应位点是否都是预期残基；
- 该位点是否位于当前筛选得到的候选孔道内壁区域；
- 该位点在孔道轴向和径向上的空间位置。

如果不先做这一步，后续 FoldX 可能会突变错误位置。
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from npstructfeat.io import get_output_file, resolve_path


# 一字母氨基酸到三字母氨基酸的映射
# 后续可以用于解析 T232K 这类突变字符串
AA1_TO_AA3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


AA3_TO_AA1 = {v: k for k, v in AA1_TO_AA3.items()}


def get_residue_features_file(config: Dict[str, Any]) -> Path:
    """
    获取 residue_features.csv 文件路径。

    默认路径规则：
        data/processed/residue_features/<PDB_ID>_residue_features.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        residue_features.csv 路径。
    """

    residue_file = get_output_file(
        config=config,
        output_dir_key="residue_feature_dir",
        suffix="residue_features.csv",
    )

    if not residue_file.exists():
        raise FileNotFoundError(
            f"残基特征表不存在: {residue_file}\n"
            "请先运行 scripts/02_extract_residue_features.py 或 scripts/run_pipeline.py"
        )

    return residue_file


def get_inner_candidate_file(config: Dict[str, Any]) -> Path:
    """
    获取 inner_candidate_residues.csv 文件路径。

    默认路径规则：
        data/processed/inner_residues/<PDB_ID>_inner_candidate_residues.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        inner_candidate_residues.csv 路径。
    """

    inner_file = get_output_file(
        config=config,
        output_dir_key="inner_residue_dir",
        suffix="inner_candidate_residues.csv",
    )

    if not inner_file.exists():
        raise FileNotFoundError(
            f"候选孔道内壁残基表不存在: {inner_file}\n"
            "请先运行 scripts/03_select_inner_residues.py 或 scripts/run_pipeline.py"
        )

    return inner_file


def get_mutation_sites_output_dir(config: Dict[str, Any]) -> Path:
    """
    获取 mutation site 检查结果输出目录。

    当前固定输出到：
        data/processed/mutation_sites/

    如果目录不存在，则自动创建。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    Path
        mutation_sites 输出目录。
    """

    output_dir = resolve_path("data/processed/mutation_sites")
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def parse_simple_mutation(mutation_id: str) -> Dict[str, Any]:
    """
    解析简单点突变字符串。

    示例
    ----
    T232K

    解析结果：
    {
        "original_aa_1": "T",
        "residue_number": 232,
        "target_aa_1": "K",
        "original_aa_3": "THR",
        "target_aa_3": "LYS"
    }

    参数
    ----
    mutation_id : str
        突变字符串，例如 "T232K"。

    返回
    ----
    dict
        解析结果。

    注意
    ----
    当前只支持最简单的单点突变格式：
        原始一字母氨基酸 + 位点编号 + 目标一字母氨基酸

    例如：
        T232K
        K238Q
        R220A

    不支持：
        T232K/K238Q
        del
        insertion
        多位点复杂写法
    """

    mutation_id = mutation_id.strip().upper()

    if len(mutation_id) < 3:
        raise ValueError(f"突变字符串过短，无法解析: {mutation_id}")

    original_aa_1 = mutation_id[0]
    target_aa_1 = mutation_id[-1]
    position_str = mutation_id[1:-1]

    if original_aa_1 not in AA1_TO_AA3:
        raise ValueError(f"未知原始氨基酸一字母缩写: {original_aa_1}")

    if target_aa_1 not in AA1_TO_AA3:
        raise ValueError(f"未知目标氨基酸一字母缩写: {target_aa_1}")

    if not position_str.isdigit():
        raise ValueError(
            f"突变位点编号无法解析: {position_str}\n"
            f"当前只支持类似 T232K 的简单格式。"
        )

    residue_number = int(position_str)

    return {
        "mutation_id": mutation_id,
        "original_aa_1": original_aa_1,
        "target_aa_1": target_aa_1,
        "original_aa_3": AA1_TO_AA3[original_aa_1],
        "target_aa_3": AA1_TO_AA3[target_aa_1],
        "residue_number": residue_number,
    }


def load_residue_and_inner_tables(
    config: Dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    读取 residue_features.csv 和 inner_candidate_residues.csv。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    residue_df, inner_df : tuple[pd.DataFrame, pd.DataFrame]
        残基特征表和候选孔道内壁残基表。
    """

    residue_file = get_residue_features_file(config)
    inner_file = get_inner_candidate_file(config)

    residue_df = pd.read_csv(residue_file)
    inner_df = pd.read_csv(inner_file)

    return residue_df, inner_df


def check_mutation_site(
    config: Dict[str, Any],
    mutation_id: str,
    chains: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    检查指定突变位点在 WT 结构中的状态。

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变字符串，例如：
        "T232K"

    chains : list[str] | None
        需要检查的链。
        如果为 None，则检查 residue_features.csv 中所有链。

    返回
    ----
    site_df : pandas.DataFrame
        突变位点检查结果。

    输出表每一行对应一条链上的突变位点。

    主要检查：
    1. 该链上 residue_number 是否存在；
    2. 实际残基是否等于预期残基；
    3. 该位点是否在候选孔道内壁残基表中；
    4. 该位点的 radial_distance、z、z_norm 等几何信息。
    """

    mutation_info = parse_simple_mutation(mutation_id)

    residue_number = mutation_info["residue_number"]
    expected_residue = mutation_info["original_aa_3"]
    target_residue = mutation_info["target_aa_3"]

    residue_df, inner_df = load_residue_and_inner_tables(config)

    pdb_id = config["input"]["pdb_id"].upper()
    nanopore_id = config["input"].get("nanopore_id", "")

    if chains is None:
        chains = sorted(residue_df["chain_id"].unique().tolist())

    rows = []

    # 为了快速判断某个 chain/residue 是否在 inner candidate 表中，
    # 构造一个集合，里面存放 (chain_id, residue_number)
    inner_key_set = set(
        zip(
            inner_df["chain_id"].astype(str),
            inner_df["residue_number"].astype(int),
        )
    )

    # inner_df 中包含 radial_distance、z_norm 等几何字段；
    # residue_df 中只有 x/y/z 和理化字段。
    # 因此如果某个位点属于 inner candidate，就优先从 inner_df 取几何信息。
    for chain_id in chains:
        chain_id = str(chain_id)

        residue_hit = residue_df[
            (residue_df["chain_id"].astype(str) == chain_id)
            & (residue_df["residue_number"].astype(int) == residue_number)
        ]

        if len(residue_hit) == 0:
            # 当前链中没有找到该编号残基
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "nanopore_id": nanopore_id,
                    "mutation_id": mutation_info["mutation_id"],
                    "chain_id": chain_id,
                    "residue_number": residue_number,
                    "expected_residue": expected_residue,
                    "observed_residue": None,
                    "target_residue": target_residue,
                    "is_site_found": False,
                    "is_expected_match": False,
                    "is_inner_candidate": False,
                    "x": None,
                    "y": None,
                    "z": None,
                    "radial_distance": None,
                    "z_norm": None,
                    "charge": None,
                    "hydrophobicity": None,
                    "is_polar": None,
                    "is_aromatic": None,
                    "note": "Residue number not found in this chain.",
                }
            )
            continue

        if len(residue_hit) > 1:
            # 正常情况下同一条链同一 residue_number 应该只有一行。
            # 如果有多个，通常可能涉及 insertion code。
            note = "Multiple residues found for this chain and residue number."
        else:
            note = ""

        residue_row = residue_hit.iloc[0]

        observed_residue = residue_row["residue_name"]
        is_expected_match = observed_residue == expected_residue

        key = (chain_id, residue_number)
        is_inner_candidate = key in inner_key_set

        # 默认使用 residue_features 中的坐标和理化性质
        x = residue_row.get("x", None)
        y = residue_row.get("y", None)
        z = residue_row.get("z", None)

        charge = residue_row.get("charge", None)
        hydrophobicity = residue_row.get("hydrophobicity", None)
        is_polar = residue_row.get("is_polar", None)
        is_aromatic = residue_row.get("is_aromatic", None)

        radial_distance = None
        z_norm = None

        # 如果该位点属于 inner candidate，则从 inner_df 获取几何增强信息
        if is_inner_candidate:
            inner_hit = inner_df[
                (inner_df["chain_id"].astype(str) == chain_id)
                & (inner_df["residue_number"].astype(int) == residue_number)
            ]

            if len(inner_hit) > 0:
                inner_row = inner_hit.iloc[0]
                radial_distance = inner_row.get("radial_distance", None)
                z_norm = inner_row.get("z_norm", None)

                # 使用 inner_df 中的 z，保持和几何表一致
                z = inner_row.get("z", z)

        rows.append(
            {
                "pdb_id": pdb_id,
                "nanopore_id": nanopore_id,
                "mutation_id": mutation_info["mutation_id"],
                "chain_id": chain_id,
                "residue_number": residue_number,
                "expected_residue": expected_residue,
                "observed_residue": observed_residue,
                "target_residue": target_residue,
                "is_site_found": True,
                "is_expected_match": bool(is_expected_match),
                "is_inner_candidate": bool(is_inner_candidate),
                "x": x,
                "y": y,
                "z": z,
                "radial_distance": radial_distance,
                "z_norm": z_norm,
                "charge": charge,
                "hydrophobicity": hydrophobicity,
                "is_polar": is_polar,
                "is_aromatic": is_aromatic,
                "note": note,
            }
        )

    site_df = pd.DataFrame(rows)

    return site_df


def save_mutation_site_check(
    config: Dict[str, Any],
    mutation_id: str,
    chains: Optional[List[str]] = None,
) -> Path:
    """
    保存突变位点检查结果。

    输出路径：
        data/processed/mutation_sites/<PDB_ID>_<MUTATION_ID>_site_check.csv

    例如：
        data/processed/mutation_sites/5JZT_T232K_site_check.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    mutation_id : str
        突变字符串，例如 "T232K"。

    chains : list[str] | None
        需要检查的链。

    返回
    ----
    Path
        输出文件路径。
    """

    site_df = check_mutation_site(
        config=config,
        mutation_id=mutation_id,
        chains=chains,
    )

    pdb_id = config["input"]["pdb_id"].upper()
    mutation_id = mutation_id.strip().upper()

    output_dir = get_mutation_sites_output_dir(config)
    output_file = output_dir / f"{pdb_id}_{mutation_id}_site_check.csv"

    site_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file


def summarize_mutation_site_check(site_df: pd.DataFrame) -> Dict[str, Any]:
    """
    对突变位点检查结果做简要统计。

    参数
    ----
    site_df : pandas.DataFrame
        check_mutation_site() 返回的结果表。

    返回
    ----
    summary : dict
        简要统计结果。
    """

    total_chains = len(site_df)

    site_found_count = int(site_df["is_site_found"].sum())
    expected_match_count = int(site_df["is_expected_match"].sum())
    inner_candidate_count = int(site_df["is_inner_candidate"].sum())

    summary = {
        "total_chains_checked": total_chains,
        "site_found_count": site_found_count,
        "expected_match_count": expected_match_count,
        "inner_candidate_count": inner_candidate_count,
        "all_sites_found": site_found_count == total_chains,
        "all_expected_match": expected_match_count == total_chains,
        "all_inner_candidate": inner_candidate_count == total_chains,
    }

    return summary