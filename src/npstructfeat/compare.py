"""
compare.py

本模块用于比较 WT 和 mutant 的结构特征差异。

当前主要用于：
1. 读取 WT nanopore-level feature table；
2. 读取 mutant nanopore-level feature table；
3. 对数值型特征计算 mutant - WT；
4. 输出 delta feature table。

为什么需要这个模块？
对于 T232K 这类突变体，单独看 mutant 特征意义有限。
更重要的是定量分析突变相对于 WT 造成了哪些结构特征变化。
"""

from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import pandas as pd

from npstructfeat.io import resolve_path
from npstructfeat.mutation import normalize_mutation_label


def read_single_row_feature_table(feature_file: str | Path) -> pd.Series:
    """
    读取单行 nanopore-level feature table。

    参数
    ----
    feature_file : str | Path
        特征表路径。

    返回
    ----
    row : pandas.Series
        第一行特征。

    注意
    ----
    当前 nanopore_structure_features.csv 默认每个结构一行。
    如果表中不是一行，说明输入可能不是当前预期格式。
    """

    feature_file = Path(feature_file)

    if not feature_file.exists():
        raise FileNotFoundError(f"特征表不存在: {feature_file}")

    df = pd.read_csv(feature_file)

    if len(df) == 0:
        raise ValueError(f"特征表为空: {feature_file}")

    if len(df) > 1:
        raise ValueError(
            f"当前函数预期特征表只有一行，但得到 {len(df)} 行: {feature_file}"
        )

    return df.iloc[0]


def is_numeric_value(value: Any) -> bool:
    """
    判断一个值是否可以视为数值。

    参数
    ----
    value : Any
        任意输入值。

    返回
    ----
    bool
        如果可以转换成 float，则返回 True。
    """

    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def compare_feature_rows(
    wt_row: pd.Series,
    mutant_row: pd.Series,
    feature_names: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    比较 WT 和 mutant 的特征行。

    参数
    ----
    wt_row : pandas.Series
        WT 特征行。

    mutant_row : pandas.Series
        mutant 特征行。

    feature_names : Iterable[str] | None
        需要比较的特征名。
        如果为 None，则自动比较两行中共同存在的数值型字段。

    返回
    ----
    delta_df : pandas.DataFrame
        差异特征表。

    输出字段
    ----
    - feature_name
    - wt_value
    - mutant_value
    - delta_value
    - relative_delta
    - comparison_note

    说明
    ----
    delta_value = mutant_value - wt_value

    relative_delta = delta_value / abs(wt_value)

    如果 wt_value 为 0，则 relative_delta 为空。
    """

    if feature_names is None:
        common_features = sorted(set(wt_row.index) & set(mutant_row.index))

        selected_features = []
        for feature in common_features:
            wt_value = wt_row[feature]
            mutant_value = mutant_row[feature]

            if is_numeric_value(wt_value) and is_numeric_value(mutant_value):
                selected_features.append(feature)
    else:
        selected_features = list(feature_names)

    rows = []

    for feature in selected_features:
        if feature not in wt_row.index:
            rows.append(
                {
                    "feature_name": feature,
                    "wt_value": None,
                    "mutant_value": None,
                    "delta_value": None,
                    "relative_delta": None,
                    "comparison_note": "Feature missing in WT row.",
                }
            )
            continue

        if feature not in mutant_row.index:
            rows.append(
                {
                    "feature_name": feature,
                    "wt_value": None,
                    "mutant_value": None,
                    "delta_value": None,
                    "relative_delta": None,
                    "comparison_note": "Feature missing in mutant row.",
                }
            )
            continue

        wt_value_raw = wt_row[feature]
        mutant_value_raw = mutant_row[feature]

        if not is_numeric_value(wt_value_raw) or not is_numeric_value(mutant_value_raw):
            rows.append(
                {
                    "feature_name": feature,
                    "wt_value": wt_value_raw,
                    "mutant_value": mutant_value_raw,
                    "delta_value": None,
                    "relative_delta": None,
                    "comparison_note": "Non-numeric feature skipped.",
                }
            )
            continue

        wt_value = float(wt_value_raw)
        mutant_value = float(mutant_value_raw)

        delta_value = mutant_value - wt_value

        if wt_value == 0:
            relative_delta = None
            note = "WT value is zero; relative_delta is undefined."
        else:
            relative_delta = delta_value / abs(wt_value)
            note = ""

        rows.append(
            {
                "feature_name": feature,
                "wt_value": wt_value,
                "mutant_value": mutant_value,
                "delta_value": delta_value,
                "relative_delta": relative_delta,
                "comparison_note": note,
            }
        )

    delta_df = pd.DataFrame(rows)

    return delta_df


def add_comparison_metadata(
    delta_df: pd.DataFrame,
    wt_id: str,
    mutant_id: str,
    template_pdb_id: str,
    mutation_id: str,
    modeling_method: str = "FoldX",
) -> pd.DataFrame:
    """
    为差异特征表添加元信息。

    参数
    ----
    delta_df : pandas.DataFrame
        差异特征表。

    wt_id : str
        WT nanopore ID，例如 aerolysin_WT。

    mutant_id : str
        mutant nanopore ID，例如 aerolysin_T232K。

    template_pdb_id : str
        模板 PDB ID，例如 5JZT。

    mutation_id : str
        突变 ID，例如 T232K。

    modeling_method : str
        突变建模方法，例如 FoldX。

    返回
    ----
    pandas.DataFrame
        添加元信息后的表。
    """

    df = delta_df.copy()

    df.insert(0, "template_pdb_id", template_pdb_id)
    df.insert(1, "wt_id", wt_id)
    df.insert(2, "mutant_id", mutant_id)
    df.insert(3, "mutation_id", mutation_id)
    df.insert(4, "modeling_method", modeling_method)

    return df


def compare_wt_mutant_feature_files(
    wt_feature_file: str | Path,
    mutant_feature_file: str | Path,
    wt_id: str,
    mutant_id: str,
    template_pdb_id: str,
    mutation_id: str,
    modeling_method: str = "FoldX",
    feature_names: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    比较 WT 和 mutant 的 nanopore-level feature files。

    参数
    ----
    wt_feature_file : str | Path
        WT nanopore_structure_features.csv 路径。

    mutant_feature_file : str | Path
        mutant nanopore_structure_features.csv 路径。

    wt_id : str
        WT ID。

    mutant_id : str
        mutant ID。

    template_pdb_id : str
        模板 PDB ID。

    mutation_id : str
        突变 ID。

    modeling_method : str
        建模方法。

    feature_names : Iterable[str] | None
        指定比较的特征名。
        如果 None，则自动比较共同的数值型特征。

    返回
    ----
    pandas.DataFrame
        差异特征表。
    """

    wt_row = read_single_row_feature_table(wt_feature_file)
    mutant_row = read_single_row_feature_table(mutant_feature_file)

    delta_df = compare_feature_rows(
        wt_row=wt_row,
        mutant_row=mutant_row,
        feature_names=feature_names,
    )

    delta_df = add_comparison_metadata(
        delta_df=delta_df,
        wt_id=wt_id,
        mutant_id=mutant_id,
        template_pdb_id=template_pdb_id,
        mutation_id=mutation_id,
        modeling_method=modeling_method,
    )

    return delta_df


def get_delta_features_output_dir() -> Path:
    """
    获取 delta features 输出目录。

    默认：
        data/processed/delta_features/

    返回
    ----
    Path
        输出目录。
    """

    output_dir = resolve_path("data/processed/delta_features")
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def save_delta_features(
    delta_df: pd.DataFrame,
    template_pdb_id: str,
    mutation_id: str,
) -> Path:
    """
    保存 WT-mutant delta features。

    输出文件：
        data/processed/delta_features/<PDB_ID>_WT_vs_<MUTATION_ID>_delta_features.csv

    示例：
        data/processed/delta_features/5JZT_WT_vs_T232K_delta_features.csv

    参数
    ----
    delta_df : pandas.DataFrame
        差异特征表。

    template_pdb_id : str
        模板 PDB ID。

    mutation_id : str
        突变 ID。

    返回
    ----
    Path
        输出文件路径。
    """

    output_dir = get_delta_features_output_dir()

    template_pdb_id = template_pdb_id.upper()
    mutation_label = normalize_mutation_label(mutation_id)

    # ?????????????? "/" ????????
    output_file = output_dir / f"{template_pdb_id}_WT_vs_{mutation_label}_delta_features.csv"

    delta_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    return output_file
