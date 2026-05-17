"""
sliding_window_pore_profile.py

沿纳米孔 z 轴方向进行滑动窗口分析。

功能：
1. 读取 candidate inner residues 表；
2. 沿 z 方向以窗口长度 L 和步长 s 滑动；
3. 统计每个窗口内的残基数量、电荷、疏水性、极性比例、芳香比例等；
4. 输出 sliding-window pore profile 表；
5. 输出滑动窗口过程示意图；
6. 输出沿 z 轴变化的 pore-profile 特征曲线图。

示例运行：
    python analysis/sliding_window_pore_profile.py

指定 T232K：
    python analysis/sliding_window_pore_profile.py ^
        --input data/processed/inner_residues/5JZT_T232K_inner_candidate_residues.csv ^
        --prefix 5JZT_T232K

指定窗口长度和步长：
    python analysis/sliding_window_pore_profile.py --window-length 10 --step-size 2
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def resolve_path(path: str | Path) -> Path:
    """
    将相对路径解析为工程根目录下的绝对路径。

    假设本脚本位于：
        analysis/sliding_window_pore_profile.py

    工程根目录为：
        analysis/ 的上一级目录
    """

    path = Path(path)

    if path.is_absolute():
        return path

    project_root = Path(__file__).resolve().parents[1]
    return project_root / path


def ensure_dir(path: str | Path) -> Path:
    """
    创建目录并返回 Path。
    """

    path = resolve_path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_inner_residues(input_file: str | Path) -> pd.DataFrame:
    """
    读取 candidate inner residues 表。

    必须包含字段：
        z
        radial_distance
        charge
        hydrophobicity
        is_polar
        is_aromatic
        residue_name
        chain_id
        residue_number
    """

    input_file = resolve_path(input_file)

    if not input_file.exists():
        raise FileNotFoundError(f"输入文件不存在: {input_file}")

    df = pd.read_csv(input_file)

    required_cols = [
        "z",
        "radial_distance",
        "charge",
        "hydrophobicity",
        "is_polar",
        "is_aromatic",
        "residue_name",
        "chain_id",
        "residue_number",
    ]

    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(
            f"输入表缺少必要字段: {missing_cols}\n"
            f"当前文件: {input_file}"
        )

    return df


def build_sliding_windows(
    z_min: float,
    z_max: float,
    window_length: float,
    step_size: float,
) -> pd.DataFrame:
    """
    根据 z_min, z_max, window_length, step_size 构建滑动窗口。

    参数
    ----
    z_min : float
        z 方向最小值。

    z_max : float
        z 方向最大值。

    window_length : float
        窗口长度，单位 Å。

    step_size : float
        步长，单位 Å。

    返回
    ----
    windows_df : pandas.DataFrame
        每一行对应一个窗口。
    """

    if window_length <= 0:
        raise ValueError("window_length 必须大于 0")

    if step_size <= 0:
        raise ValueError("step_size 必须大于 0")

    rows = []

    start = z_min
    window_id = 0

    while start <= z_max:
        end = start + window_length

        rows.append(
            {
                "window_id": window_id,
                "z_start": start,
                "z_end": end,
                "z_center": (start + end) / 2.0,
            }
        )

        start += step_size
        window_id += 1

        # 如果窗口起点已经超过 z_max，就停止
        if start > z_max:
            break

    return pd.DataFrame(rows)


def summarize_window(
    window_df: pd.DataFrame,
    z_start: float,
    z_end: float,
) -> dict:
    """
    对单个窗口内的残基进行统计。

    窗口定义：
        z_start <= z < z_end

    返回一个 dict，后续作为一行输出。
    """

    residue_count = len(window_df)

    if residue_count == 0:
        return {
            "residue_count": 0,
            "positive_count": 0,
            "negative_count": 0,
            "neutral_count": 0,
            "net_charge": 0,
            "mean_charge": None,
            "mean_hydrophobicity": None,
            "std_hydrophobicity": None,
            "polar_count": 0,
            "polar_ratio": None,
            "aromatic_count": 0,
            "aromatic_ratio": None,
            "mean_radial_distance": None,
            "min_radial_distance": None,
            "max_radial_distance": None,
            "residue_type_count": 0,
            "chain_count": 0,
        }

    charge = window_df["charge"].astype(float)

    positive_count = int((charge > 0).sum())
    negative_count = int((charge < 0).sum())
    neutral_count = int((charge == 0).sum())

    polar_count = int(window_df["is_polar"].astype(int).sum())
    aromatic_count = int(window_df["is_aromatic"].astype(int).sum())

    return {
        "residue_count": int(residue_count),
        "positive_count": positive_count,
        "negative_count": negative_count,
        "neutral_count": neutral_count,
        "net_charge": float(charge.sum()),
        "mean_charge": float(charge.mean()),
        "mean_hydrophobicity": float(window_df["hydrophobicity"].mean()),
        "std_hydrophobicity": float(window_df["hydrophobicity"].std(ddof=0)),
        "polar_count": polar_count,
        "polar_ratio": polar_count / residue_count,
        "aromatic_count": aromatic_count,
        "aromatic_ratio": aromatic_count / residue_count,
        "mean_radial_distance": float(window_df["radial_distance"].mean()),
        "min_radial_distance": float(window_df["radial_distance"].min()),
        "max_radial_distance": float(window_df["radial_distance"].max()),
        "residue_type_count": int(window_df["residue_name"].nunique()),
        "chain_count": int(window_df["chain_id"].nunique()),
    }


def compute_sliding_window_profile(
    inner_df: pd.DataFrame,
    window_length: float,
    step_size: float,
) -> pd.DataFrame:
    """
    计算 sliding-window pore profile。

    参数
    ----
    inner_df : pandas.DataFrame
        candidate inner residues 表。

    window_length : float
        窗口长度，单位 Å。

    step_size : float
        步长，单位 Å。

    返回
    ----
    profile_df : pandas.DataFrame
        每一行对应一个 z-window 的统计特征。
    """

    z_min = float(inner_df["z"].min())
    z_max = float(inner_df["z"].max())

    windows_df = build_sliding_windows(
        z_min=z_min,
        z_max=z_max,
        window_length=window_length,
        step_size=step_size,
    )

    rows = []

    for _, win in windows_df.iterrows():
        z_start = float(win["z_start"])
        z_end = float(win["z_end"])

        window_residues = inner_df[
            (inner_df["z"] >= z_start)
            & (inner_df["z"] < z_end)
        ].copy()

        summary = summarize_window(
            window_df=window_residues,
            z_start=z_start,
            z_end=z_end,
        )

        row = {
            "window_id": int(win["window_id"]),
            "z_start": z_start,
            "z_end": z_end,
            "z_center": float(win["z_center"]),
            "window_length": window_length,
            "step_size": step_size,
        }

        row.update(summary)
        rows.append(row)

    profile_df = pd.DataFrame(rows)

    return profile_df


def plot_sliding_window_process(
    inner_df: pd.DataFrame,
    profile_df: pd.DataFrame,
    output_file: str | Path,
    max_windows_to_draw: int = 12,
) -> None:
    """
    绘制滑动窗口过程示意图。

    图中：
    - 蓝色散点表示 candidate inner residues 的 z-radial 分布；
    - 半透明矩形表示沿 z 轴滑动的窗口；
    - 为避免太密，只绘制前 max_windows_to_draw 个窗口。
    """

    output_file = Path(output_file)

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.scatter(
        inner_df["z"],
        inner_df["radial_distance"],
        s=16,
        alpha=0.35,
        label="Candidate inner residues",
    )

    y_min = float(inner_df["radial_distance"].min()) - 0.5
    y_max = float(inner_df["radial_distance"].max()) + 0.5

    draw_df = profile_df.head(max_windows_to_draw)

    for _, row in draw_df.iterrows():
        z_start = float(row["z_start"])
        z_end = float(row["z_end"])

        ax.axvspan(
            z_start,
            z_end,
            alpha=0.12,
        )

        ax.text(
            (z_start + z_end) / 2.0,
            y_max,
            f"W{int(row['window_id'])}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_xlabel("z coordinate (Å)")
    ax.set_ylabel("Radial distance to pore axis (Å)")
    ax.set_title("Sliding-window process along the pore axis")
    ax.set_ylim(y_min, y_max + 1.0)
    ax.grid(linestyle="--", alpha=0.35)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def plot_profile_features(
    profile_df: pd.DataFrame,
    output_file: str | Path,
) -> None:
    """
    绘制 sliding-window 特征曲线图。

    为避免多轴混乱，这里使用 3 张独立图拼接在一个 figure 中：
    1. residue_count
    2. net_charge
    3. mean_hydrophobicity, polar_ratio, aromatic_ratio

    注意：
    这里使用 subplots 是为了在同一张输出图中展示 pore profile 的多维趋势。
    """

    output_file = Path(output_file)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    x = profile_df["z_center"]

    axes[0].plot(x, profile_df["residue_count"], marker="o")
    axes[0].set_ylabel("Residue count")
    axes[0].set_title("Sliding-window pore profile")
    axes[0].grid(linestyle="--", alpha=0.35)

    axes[1].plot(x, profile_df["net_charge"], marker="o")
    axes[1].axhline(0, linewidth=1)
    axes[1].set_ylabel("Net charge")
    axes[1].grid(linestyle="--", alpha=0.35)

    axes[2].plot(x, profile_df["mean_hydrophobicity"], marker="o", label="Mean hydrophobicity")
    axes[2].plot(x, profile_df["polar_ratio"], marker="o", label="Polar ratio")
    axes[2].plot(x, profile_df["aromatic_ratio"], marker="o", label="Aromatic ratio")
    axes[2].set_xlabel("Window center z coordinate (Å)")
    axes[2].set_ylabel("Feature value")
    axes[2].grid(linestyle="--", alpha=0.35)
    axes[2].legend()

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    """
    解析命令行参数。
    """

    parser = argparse.ArgumentParser(
        description="Build sliding-window pore profile along the nanopore z-axis."
    )

    parser.add_argument(
        "--input",
        type=str,
        default="data/processed/inner_residues/5JZT_inner_candidate_residues.csv",
        help="Input candidate inner residue CSV file.",
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default="5JZT",
        help="Output prefix, e.g., 5JZT or 5JZT_T232K.",
    )

    parser.add_argument(
        "--window-length",
        type=float,
        default=10.0,
        help="Sliding window length along z-axis, in Å. Default: 10.0",
    )

    parser.add_argument(
        "--step-size",
        type=float,
        default=2.0,
        help="Sliding step size along z-axis, in Å. Default: 2.0",
    )

    parser.add_argument(
        "--profile-output-dir",
        type=str,
        default="data/processed/pore_profiles",
        help="Output directory for sliding-window profile CSV.",
    )

    parser.add_argument(
        "--figure-output-dir",
        type=str,
        default="outputs/figures",
        help="Output directory for figures.",
    )

    return parser.parse_args()


def main() -> None:
    """
    主函数。
    """

    args = parse_args()

    input_file = resolve_path(args.input)
    profile_output_dir = ensure_dir(args.profile_output_dir)
    figure_output_dir = ensure_dir(args.figure_output_dir)

    print("Sliding-window pore profile analysis")
    print("=" * 60)
    print(f"Input file:      {input_file}")
    print(f"Prefix:          {args.prefix}")
    print(f"Window length:   {args.window_length} Å")
    print(f"Step size:       {args.step_size} Å")
    print(f"Profile out dir: {profile_output_dir}")
    print(f"Figure out dir:  {figure_output_dir}")

    inner_df = load_inner_residues(input_file)

    print("\nInput summary:")
    print(f"Candidate residue count: {len(inner_df)}")
    print(f"z_min: {inner_df['z'].min():.3f}")
    print(f"z_max: {inner_df['z'].max():.3f}")
    print(f"z_range: {inner_df['z'].max() - inner_df['z'].min():.3f} Å")

    profile_df = compute_sliding_window_profile(
        inner_df=inner_df,
        window_length=args.window_length,
        step_size=args.step_size,
    )

    profile_file = profile_output_dir / (
        f"{args.prefix}_sliding_window_profile_"
        f"L{args.window_length:g}_S{args.step_size:g}.csv"
    )

    profile_df.to_csv(profile_file, index=False, encoding="utf-8-sig")

    process_fig = figure_output_dir / (
        f"{args.prefix}_sliding_window_process_"
        f"L{args.window_length:g}_S{args.step_size:g}.png"
    )

    profile_fig = figure_output_dir / (
        f"{args.prefix}_sliding_window_profile_features_"
        f"L{args.window_length:g}_S{args.step_size:g}.png"
    )

    plot_sliding_window_process(
        inner_df=inner_df,
        profile_df=profile_df,
        output_file=process_fig,
    )

    plot_profile_features(
        profile_df=profile_df,
        output_file=profile_fig,
    )

    print("\nSliding-window profile preview:")
    display_cols = [
        "window_id",
        "z_start",
        "z_end",
        "z_center",
        "residue_count",
        "net_charge",
        "mean_hydrophobicity",
        "polar_ratio",
        "aromatic_ratio",
        "mean_radial_distance",
    ]

    print(
        profile_df[display_cols].head(10).to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print("\n输出文件：")
    print(profile_file)
    print(process_fig)
    print(profile_fig)

    print("\n解释建议：")
    print("1. sliding_window_profile.csv 可作为 pore-profile 结构特征表。")
    print("2. sliding_window_process.png 展示窗口如何沿 z 轴滑动。")
    print("3. sliding_window_profile_features.png 展示残基数量、电荷、疏水性等沿 z 轴的变化。")
    print("4. 如果要比较 WT 和 T232K，应分别运行该脚本，并比较两个 profile CSV。")


if __name__ == "__main__":
    main()