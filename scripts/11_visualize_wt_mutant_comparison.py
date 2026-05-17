"""
11_visualize_wt_mutant_comparison.py

本脚本用于可视化比较 WT aerolysin 与 FoldX-modeled T232K mutant 的结构特征差异。

当前默认输入：
1. WT nanopore-level feature table
   data/processed/nanopore_features/5JZT_nanopore_structure_features.csv

2. T232K nanopore-level feature table
   data/processed/nanopore_features/5JZT_T232K_nanopore_structure_features.csv

3. WT-mutant delta feature table
   data/processed/delta_features/5JZT_WT_vs_T232K_delta_features.csv

4. WT inner candidate residues
   data/processed/inner_residues/5JZT_inner_candidate_residues.csv

5. T232K mutation site check table
   data/processed/mutation_sites/5JZT_T232K_site_check.csv

输出图像：
1. WT vs T232K key features
2. WT-mutant delta features
3. WT vs T232K charge composition
4. T232 site position in z-radial space

运行方式：
    python scripts/11_visualize_wt_mutant_comparison.py

如果找不到 npstructfeat：
    $env:PYTHONPATH="src"
    python scripts/11_visualize_wt_mutant_comparison.py
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from npstructfeat.io import resolve_path


def read_single_row_csv(file_path: str | Path) -> pd.Series:
    """
    读取单行 CSV 表格。

    当前 nanopore_structure_features.csv 默认每个结构一行。
    """

    file_path = resolve_path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"文件不存在: {file_path}")

    df = pd.read_csv(file_path)

    if len(df) != 1:
        raise ValueError(
            f"预期文件只有一行，但实际有 {len(df)} 行: {file_path}"
        )

    return df.iloc[0]


def ensure_figure_dir(output_dir: str | Path) -> Path:
    """
    创建图像输出目录。
    """

    output_dir = resolve_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def plot_key_features(
    wt_row: pd.Series,
    mutant_row: pd.Series,
    output_file: str | Path,
) -> None:
    """
    绘制 WT vs T232K 关键结构特征柱状图。

    该图用于直观比较：
    - 内壁净电荷
    - 正电残基数量
    - 负电残基数量
    - 平均疏水性
    - 极性比例
    - 芳香比例

    注意：
    不同特征数值尺度不同，因此该图主要用于定性比较。
    """

    features = [
        "inner_net_charge",
        "inner_positive_residue_count",
        "inner_negative_residue_count",
        "inner_mean_hydrophobicity",
        "inner_polar_ratio",
        "inner_aromatic_ratio",
    ]

    labels = [
        "Net charge",
        "Positive count",
        "Negative count",
        "Mean hydrophobicity",
        "Polar ratio",
        "Aromatic ratio",
    ]

    wt_values = [float(wt_row[f]) for f in features]
    mutant_values = [float(mutant_row[f]) for f in features]

    x = range(len(features))
    width = 0.36

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.bar([i - width / 2 for i in x], wt_values, width, label="WT")
    ax.bar([i + width / 2 for i in x], mutant_values, width, label="T232K")

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Feature value")
    ax.set_title("WT vs T232K: Key nanopore structural features")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_delta_features(
    delta_df: pd.DataFrame,
    output_file: str | Path,
) -> None:
    """
    绘制 WT-mutant delta feature 条形图。

    delta_value = mutant_value - wt_value

    该图用于突出 T232K 相对于 WT 的变化方向和变化幅度。
    """

    key_features = [
        "inner_net_charge",
        "inner_positive_residue_count",
        "inner_negative_residue_count",
        "inner_neutral_residue_count",
        "inner_mean_charge",
        "inner_mean_hydrophobicity",
        "inner_polar_ratio",
        "inner_aromatic_ratio",
        "inner_candidate_residue_count",
        "radial_distance_mean_A",
        "pore_region_length_approx_A",
    ]

    feature_labels = {
        "inner_net_charge": "Net charge",
        "inner_positive_residue_count": "Positive count",
        "inner_negative_residue_count": "Negative count",
        "inner_neutral_residue_count": "Neutral count",
        "inner_mean_charge": "Mean charge",
        "inner_mean_hydrophobicity": "Mean hydrophobicity",
        "inner_polar_ratio": "Polar ratio",
        "inner_aromatic_ratio": "Aromatic ratio",
        "inner_candidate_residue_count": "Candidate count",
        "radial_distance_mean_A": "Mean radial distance",
        "pore_region_length_approx_A": "Pore region length",
    }

    plot_df = delta_df[delta_df["feature_name"].isin(key_features)].copy()

    plot_df["feature_order"] = plot_df["feature_name"].apply(
        lambda x: key_features.index(x)
    )
    plot_df = plot_df.sort_values("feature_order")

    labels = [feature_labels[f] for f in plot_df["feature_name"]]
    values = plot_df["delta_value"].astype(float).tolist()

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.bar(labels, values)
    ax.axhline(0, linewidth=1)

    ax.set_ylabel("Delta value (T232K - WT)")
    ax.set_title("WT vs T232K: Delta structural features")
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_charge_composition(
    wt_row: pd.Series,
    mutant_row: pd.Series,
    output_file: str | Path,
) -> None:
    """
    绘制 WT 与 T232K 的候选孔道区域电荷组成。

    重点展示：
    - positive residue count
    - negative residue count
    - neutral residue count

    对 T232K，预期结果是：
    - positive count +7
    - neutral count -7
    - negative count unchanged
    """

    groups = ["WT", "T232K"]

    positive = [
        float(wt_row["inner_positive_residue_count"]),
        float(mutant_row["inner_positive_residue_count"]),
    ]

    negative = [
        float(wt_row["inner_negative_residue_count"]),
        float(mutant_row["inner_negative_residue_count"]),
    ]

    neutral = [
        float(wt_row["inner_neutral_residue_count"]),
        float(mutant_row["inner_neutral_residue_count"]),
    ]

    x = range(len(groups))

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(x, positive, label="Positive residues")
    ax.bar(x, negative, bottom=positive, label="Negative residues")

    bottom_neutral = [p + n for p, n in zip(positive, negative)]
    ax.bar(x, neutral, bottom=bottom_neutral, label="Neutral residues")

    ax.set_xticks(list(x))
    ax.set_xticklabels(groups)
    ax.set_ylabel("Residue count")
    ax.set_title("Charge composition of candidate pore-region residues")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_t232_z_radial_position(
    inner_residue_file: str | Path,
    mutation_site_file: str | Path,
    output_file: str | Path,
) -> None:
    """
    在 z-radial 分布图中高亮 T232 位点。

    背景点：
        WT candidate inner residues

    高亮点：
        A-G 七条链的 THR232 位点

    该图用于说明 T232 位于当前几何定义下的候选孔道区域内。
    """

    inner_residue_file = resolve_path(inner_residue_file)
    mutation_site_file = resolve_path(mutation_site_file)

    if not inner_residue_file.exists():
        raise FileNotFoundError(f"inner residue file 不存在: {inner_residue_file}")

    if not mutation_site_file.exists():
        raise FileNotFoundError(f"mutation site file 不存在: {mutation_site_file}")

    inner_df = pd.read_csv(inner_residue_file)
    site_df = pd.read_csv(mutation_site_file)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        inner_df["z"],
        inner_df["radial_distance"],
        s=18,
        alpha=0.35,
        label="Candidate inner residues",
    )

    ax.scatter(
        site_df["z"],
        site_df["radial_distance"],
        s=80,
        marker="*",
        label="THR232 sites",
    )

    for _, row in site_df.iterrows():
        ax.text(
            row["z"],
            row["radial_distance"] + 0.25,
            str(row["chain_id"]),
            fontsize=8,
            ha="center",
        )

    ax.axhline(20.0, linestyle="--", linewidth=1, label="20 Å threshold")

    ax.set_xlabel("z coordinate (Å)")
    ax.set_ylabel("Radial distance to pore axis (Å)")
    ax.set_title("Position of T232 sites in candidate pore-region space")
    ax.legend()
    ax.grid(linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def parse_args():
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Visualize WT vs T232K nanopore structural feature comparison."
    )

    parser.add_argument(
        "--wt-feature",
        type=str,
        default="data/processed/nanopore_features/5JZT_nanopore_structure_features.csv",
        help="WT nanopore structure feature table.",
    )

    parser.add_argument(
        "--mutant-feature",
        type=str,
        default="data/processed/nanopore_features/5JZT_T232K_nanopore_structure_features.csv",
        help="T232K nanopore structure feature table.",
    )

    parser.add_argument(
        "--delta-feature",
        type=str,
        default="data/processed/delta_features/5JZT_WT_vs_T232K_delta_features.csv",
        help="WT vs T232K delta feature table.",
    )

    parser.add_argument(
        "--inner-residue",
        type=str,
        default="data/processed/inner_residues/5JZT_inner_candidate_residues.csv",
        help="WT candidate inner residue table.",
    )

    parser.add_argument(
        "--mutation-site",
        type=str,
        default="data/processed/mutation_sites/5JZT_T232K_site_check.csv",
        help="T232K mutation site check table.",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/figures",
        help="Output figure directory.",
    )

    return parser.parse_args()


def main():
    """
    主函数。
    """

    args = parse_args()

    output_dir = ensure_figure_dir(args.output_dir)

    wt_feature_file = resolve_path(args.wt_feature)
    mutant_feature_file = resolve_path(args.mutant_feature)
    delta_feature_file = resolve_path(args.delta_feature)

    print("读取输入文件：")
    print(f"WT feature:      {wt_feature_file}")
    print(f"Mutant feature:  {mutant_feature_file}")
    print(f"Delta feature:   {delta_feature_file}")
    print(f"Inner residues:  {resolve_path(args.inner_residue)}")
    print(f"Mutation sites:  {resolve_path(args.mutation_site)}")

    wt_row = read_single_row_csv(wt_feature_file)
    mutant_row = read_single_row_csv(mutant_feature_file)

    if not delta_feature_file.exists():
        raise FileNotFoundError(
            f"delta feature file 不存在: {delta_feature_file}\n"
            "请先运行 scripts/10_compare_wt_mutant_features.py"
        )

    delta_df = pd.read_csv(delta_feature_file)

    fig1 = output_dir / "5JZT_WT_vs_T232K_key_features.png"
    fig2 = output_dir / "5JZT_WT_vs_T232K_delta_features.png"
    fig3 = output_dir / "5JZT_WT_vs_T232K_charge_composition.png"
    fig4 = output_dir / "5JZT_T232_site_z_radial_position.png"

    plot_key_features(
        wt_row=wt_row,
        mutant_row=mutant_row,
        output_file=fig1,
    )

    plot_delta_features(
        delta_df=delta_df,
        output_file=fig2,
    )

    plot_charge_composition(
        wt_row=wt_row,
        mutant_row=mutant_row,
        output_file=fig3,
    )

    plot_t232_z_radial_position(
        inner_residue_file=args.inner_residue,
        mutation_site_file=args.mutation_site,
        output_file=fig4,
    )

    print("\n可视化图像已保存：")
    print(fig1)
    print(fig2)
    print(fig3)
    print(fig4)

    print("\n图像解释：")
    print("1. key_features：比较 WT 与 T232K 的关键结构特征。")
    print("2. delta_features：展示 T232K - WT 的特征变化。")
    print("3. charge_composition：展示候选孔道区域正/负/中性残基组成变化。")
    print("4. T232_site_z_radial_position：展示 T232 位点在 z-radial 空间中的位置。")


if __name__ == "__main__":
    main()