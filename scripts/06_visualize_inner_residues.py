"""
06_visualize_inner_residues.py

本脚本用于可视化检查候选孔道内壁残基筛选结果。

为什么需要可视化？
前面的 03_select_inner_residues.py 已经基于 radial_distance <= threshold
筛选出了候选孔道内壁残基。但是仅凭 CSV 表格不够直观。

本脚本主要检查三个问题：

1. radial_distance 的整体分布是否合理？
   - 当前 20 Å 阈值切在什么位置？
   - 是否明显过宽或过窄？

2. 候选残基是否沿 z 轴连续覆盖孔道区域？
   - 如果只集中在局部 z 区域，说明筛选可能不合理；
   - 如果沿 z 轴覆盖较长区域，说明更符合孔道内壁候选区域。

3. 阈值敏感性是否平滑？
   - threshold 从 15 到 25 Å 时，候选残基数量、电荷、疏水性、极性比例是否突变？

输入文件：
    data/processed/residue_features/5JZT_residue_features.csv
    data/processed/inner_residues/5JZT_inner_candidate_residues.csv
    data/processed/nanopore_features/5JZT_threshold_sensitivity.csv

输出图像：
    outputs/figures/5JZT_radial_distance_hist.png
    outputs/figures/5JZT_z_radial_scatter.png
    outputs/figures/5JZT_threshold_sensitivity_count.png
    outputs/figures/5JZT_threshold_sensitivity_properties.png

运行方式：
    python scripts/06_visualize_inner_residues.py

或者指定配置文件：
    python scripts/06_visualize_inner_residues.py --config config/default.yaml
"""

import argparse
from pathlib import Path
from typing import Dict, Any

import matplotlib.pyplot as plt
import pandas as pd

from npstructfeat.geometry import add_basic_geometry_features
from npstructfeat.io import get_output_file, load_config, resolve_path
from npstructfeat.pore import load_residue_features
from npstructfeat.utils import print_config_summary, require_config_keys


def parse_args():
    """
    解析命令行参数。

    当前只设置一个参数：
    --config

    默认读取：
        config/default.yaml
    """

    parser = argparse.ArgumentParser(
        description="Visualize candidate pore-lining residues and threshold sensitivity."
    )

    parser.add_argument(
        "--config",
        type=str,
        default="config/default.yaml",
        help="Path to YAML config file. Default: config/default.yaml",
    )

    return parser.parse_args()


def get_figure_dir() -> Path:
    """
    获取图像输出目录。

    当前固定输出到：
        outputs/figures/

    如果目录不存在，则自动创建。

    返回
    ----
    fig_dir : Path
        图像输出目录的绝对路径。
    """

    fig_dir = resolve_path("outputs/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    return fig_dir


def get_inner_candidate_file(config: Dict[str, Any]) -> Path:
    """
    获取候选孔道内壁残基表路径。

    默认文件：
        data/processed/inner_residues/<PDB_ID>_inner_candidate_residues.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    inner_file : Path
        候选孔道内壁残基 CSV 路径。
    """

    inner_file = get_output_file(
        config=config,
        output_dir_key="inner_residue_dir",
        suffix="inner_candidate_residues.csv",
    )

    if not inner_file.exists():
        raise FileNotFoundError(
            f"候选孔道内壁残基表不存在: {inner_file}\n"
            "请先运行 scripts/03_select_inner_residues.py"
        )

    return inner_file


def get_threshold_sensitivity_file(config: Dict[str, Any]) -> Path:
    """
    获取阈值敏感性分析结果文件路径。

    默认文件：
        data/processed/nanopore_features/<PDB_ID>_threshold_sensitivity.csv

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    sensitivity_file : Path
        阈值敏感性分析 CSV 文件路径。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    sensitivity_file = resolve_path(
        f"data/processed/nanopore_features/{pdb_id}_threshold_sensitivity.csv"
    )

    if not sensitivity_file.exists():
        raise FileNotFoundError(
            f"阈值敏感性分析文件不存在: {sensitivity_file}\n"
            "请先运行 scripts/05_threshold_sensitivity.py"
        )

    return sensitivity_file


def load_inner_candidates(config: Dict[str, Any]) -> pd.DataFrame:
    """
    读取候选孔道内壁残基表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    inner_df : pandas.DataFrame
        候选孔道内壁残基表。
    """

    inner_file = get_inner_candidate_file(config)
    inner_df = pd.read_csv(inner_file)

    return inner_df


def load_threshold_sensitivity(config: Dict[str, Any]) -> pd.DataFrame:
    """
    读取阈值敏感性分析结果表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    sensitivity_df : pandas.DataFrame
        阈值敏感性分析表。
    """

    sensitivity_file = get_threshold_sensitivity_file(config)
    sensitivity_df = pd.read_csv(sensitivity_file)

    return sensitivity_df


def prepare_residue_geometry(config: Dict[str, Any]) -> pd.DataFrame:
    """
    读取 residue_features.csv，并添加几何特征。

    为什么这里重新添加几何特征？
    residue_features.csv 中只有 x/y/z 和残基理化性质，
    不一定包含 radial_distance。

    这里重新调用 add_basic_geometry_features，
    得到包含 radial_distance、theta、z_norm 的完整表。

    参数
    ----
    config : dict
        YAML 配置字典。

    返回
    ----
    residue_geo_df : pandas.DataFrame
        添加几何特征后的残基表。
    """

    residue_df = load_residue_features(config)

    center_mode = config.get("pore", {}).get("center_mode", "xy_mean")

    residue_geo_df = add_basic_geometry_features(
        residue_df=residue_df,
        center_mode=center_mode,
    )

    return residue_geo_df


def plot_radial_distance_hist(
    residue_geo_df: pd.DataFrame,
    config: Dict[str, Any],
    fig_dir: Path,
) -> Path:
    """
    绘制所有残基 radial_distance 分布直方图。

    图像目的
    ----
    观察当前 inner_radius_threshold 是否切在合理位置。

    x 轴：
        radial_distance，单位 Å

    y 轴：
        residue count

    图中会用竖线标出当前阈值，例如 20 Å。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    threshold = float(config.get("pore", {}).get("inner_radius_threshold", 20.0))

    output_file = fig_dir / f"{pdb_id}_radial_distance_hist.png"

    plt.figure(figsize=(8, 5))

    plt.hist(
        residue_geo_df["radial_distance"],
        bins=40,
        edgecolor="black",
        alpha=0.75,
    )

    # 当前阈值线
    plt.axvline(
        threshold,
        linestyle="--",
        linewidth=2,
        label=f"threshold = {threshold:.1f} Å",
    )

    plt.xlabel("Radial distance to pore axis (Å)", fontsize=12)
    plt.ylabel("Residue count", fontsize=12)
    plt.title(f"{pdb_id} radial distance distribution", fontsize=14)
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file


def plot_z_radial_scatter(
    residue_geo_df: pd.DataFrame,
    inner_df: pd.DataFrame,
    config: Dict[str, Any],
    fig_dir: Path,
) -> Path:
    """
    绘制 z - radial_distance 散点图。

    图像目的
    ----
    观察候选孔道内壁残基是否沿 z 轴连续分布。

    x 轴：
        z coordinate，单位 Å

    y 轴：
        radial_distance，单位 Å

    图像中：
    - 所有残基用较浅的点表示；
    - 候选孔道内壁残基用较深的点表示；
    - 当前 threshold 用水平虚线表示。

    注意
    ----
    这张图不是三维结构图。
    它是二维几何投影检查，用于快速判断筛选是否合理。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    threshold = float(config.get("pore", {}).get("inner_radius_threshold", 20.0))

    output_file = fig_dir / f"{pdb_id}_z_radial_scatter.png"

    plt.figure(figsize=(9, 5))

    # 所有残基
    plt.scatter(
        residue_geo_df["z"],
        residue_geo_df["radial_distance"],
        s=8,
        alpha=0.25,
        label="All residues",
    )

    # 候选孔道内壁残基
    plt.scatter(
        inner_df["z"],
        inner_df["radial_distance"],
        s=14,
        alpha=0.8,
        label="Candidate inner residues",
    )

    # 当前阈值水平线
    plt.axhline(
        threshold,
        linestyle="--",
        linewidth=2,
        label=f"threshold = {threshold:.1f} Å",
    )

    plt.xlabel("z coordinate (Å)", fontsize=12)
    plt.ylabel("Radial distance to pore axis (Å)", fontsize=12)
    plt.title(f"{pdb_id} z-radial distribution", fontsize=14)
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file


def plot_threshold_sensitivity_count(
    sensitivity_df: pd.DataFrame,
    config: Dict[str, Any],
    fig_dir: Path,
) -> Path:
    """
    绘制 threshold 与候选残基数量/比例的关系。

    图像目的
    ----
    判断候选残基数量是否随阈值平滑变化。
    如果某个阈值附近出现突变，说明该阈值附近可能纳入了大量外侧区域残基。

    当前图只画候选残基数量。
    候选比例可以从 CSV 表中查看。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    output_file = fig_dir / f"{pdb_id}_threshold_sensitivity_count.png"

    plt.figure(figsize=(8, 5))

    plt.plot(
        sensitivity_df["threshold_A"],
        sensitivity_df["inner_candidate_residue_count"],
        marker="o",
        linewidth=2,
    )

    plt.xlabel("Threshold (Å)", fontsize=12)
    plt.ylabel("Candidate residue count", fontsize=12)
    plt.title(f"{pdb_id} threshold sensitivity: residue count", fontsize=14)
    plt.tight_layout()

    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file


def plot_threshold_sensitivity_properties(
    sensitivity_df: pd.DataFrame,
    config: Dict[str, Any],
    fig_dir: Path,
) -> Path:
    """
    绘制 threshold 与关键理化性质的关系。

    图像目的
    ----
    观察不同 threshold 下：
    - inner_mean_hydrophobicity
    - inner_aromatic_ratio
    - inner_polar_ratio
    - inner_net_charge

    是否发生明显变化。

    注意
    ----
    这些指标数值尺度不同。
    为了保持图简单，当前只画 3 个相近尺度的性质：
    - 平均疏水性；
    - 芳香族比例；
    - 极性比例。

    inner_net_charge 建议从 CSV 表中直接看。
    如果把净电荷也画在同一张图上，会因为尺度不同影响阅读。
    """

    pdb_id = config["input"]["pdb_id"].upper()
    output_file = fig_dir / f"{pdb_id}_threshold_sensitivity_properties.png"

    plt.figure(figsize=(8, 5))

    plt.plot(
        sensitivity_df["threshold_A"],
        sensitivity_df["inner_mean_hydrophobicity"],
        marker="o",
        linewidth=2,
        label="Mean hydrophobicity",
    )

    plt.plot(
        sensitivity_df["threshold_A"],
        sensitivity_df["inner_aromatic_ratio"],
        marker="o",
        linewidth=2,
        label="Aromatic ratio",
    )

    plt.plot(
        sensitivity_df["threshold_A"],
        sensitivity_df["inner_polar_ratio"],
        marker="o",
        linewidth=2,
        label="Polar ratio",
    )

    plt.xlabel("Threshold (Å)", fontsize=12)
    plt.ylabel("Feature value", fontsize=12)
    plt.title(f"{pdb_id} threshold sensitivity: physicochemical features", fontsize=14)
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file


def main():
    """
    主函数。

    执行流程：
    1. 读取配置文件；
    2. 读取 residue_features.csv 并添加几何特征；
    3. 读取 candidate inner residues；
    4. 读取 threshold sensitivity 表；
    5. 输出 4 张可视化图。
    """

    args = parse_args()

    config = load_config(args.config)

    require_config_keys(
        config,
        required_keys=["project", "input", "structure", "pore", "output", "save"],
    )

    print_config_summary(config)

    fig_dir = get_figure_dir()
    print(f"\n图像输出目录: {fig_dir}")

    # 1. 准备数据
    residue_geo_df = prepare_residue_geometry(config)
    inner_df = load_inner_candidates(config)
    sensitivity_df = load_threshold_sensitivity(config)

    print("\n数据读取完成：")
    print(f"全部残基数: {len(residue_geo_df)}")
    print(f"候选孔道内壁残基数: {len(inner_df)}")
    print(f"阈值敏感性测试数量: {len(sensitivity_df)}")

    # 2. 绘图
    output_files = []

    output_files.append(
        plot_radial_distance_hist(
            residue_geo_df=residue_geo_df,
            config=config,
            fig_dir=fig_dir,
        )
    )

    output_files.append(
        plot_z_radial_scatter(
            residue_geo_df=residue_geo_df,
            inner_df=inner_df,
            config=config,
            fig_dir=fig_dir,
        )
    )

    output_files.append(
        plot_threshold_sensitivity_count(
            sensitivity_df=sensitivity_df,
            config=config,
            fig_dir=fig_dir,
        )
    )

    output_files.append(
        plot_threshold_sensitivity_properties(
            sensitivity_df=sensitivity_df,
            config=config,
            fig_dir=fig_dir,
        )
    )

    print("\n可视化图像已保存：")
    for file in output_files:
        print(file)

    print(
        "\n检查建议：\n"
        "1. radial_distance_hist：观察 20 Å 阈值是否位于合理截断位置；\n"
        "2. z_radial_scatter：观察候选残基是否沿 z 方向连续覆盖孔道区域；\n"
        "3. threshold_sensitivity_count：观察候选数量是否随阈值平滑增加；\n"
        "4. threshold_sensitivity_properties：观察理化性质是否在 20 Å 附近稳定。"
    )


if __name__ == "__main__":
    main()