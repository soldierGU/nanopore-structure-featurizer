"""
geometry.py

本文件存放与几何计算有关的函数。

当前阶段主要用于：
1. 估计纳米孔中心轴；
2. 计算残基到中心轴的径向距离；
3. 计算残基在 x-y 平面上的角度；
4. 为后续候选孔道内壁残基筛选提供几何特征。

为什么要单独写 geometry.py？
因为几何计算以后会反复用到，例如：
- 计算孔道内壁残基；
- 构建残基接触图；
- 计算孔道轴向剖面；
- 计算残基空间邻接关系。
"""

from typing import Tuple

import numpy as np
import pandas as pd


def estimate_xy_center(
    residue_df: pd.DataFrame,
    center_mode: str = "xy_mean",
) -> Tuple[float, float]:
    """
    估计孔道中心轴在 x-y 平面上的位置。

    参数
    ----
    residue_df : pandas.DataFrame
        残基级特征表。
        至少需要包含两列：
        - x
        - y

    center_mode : str
        中心估计方式。

        当前支持：
        1. "xy_mean"
           使用所有残基 Cα 坐标的 x/y 均值作为中心。

        2. "xy_median"
           使用所有残基 Cα 坐标的 x/y 中位数作为中心。

    返回
    ----
    center_x, center_y : tuple[float, float]
        估计得到的中心轴位置。

    说明
    ----
    对于高度对称的七聚体孔道，x/y 均值通常可以作为第一版近似中心。
    但如果结构没有对齐，或者孔道轴不沿 z 轴，这种方法会有偏差。
    """

    if "x" not in residue_df.columns or "y" not in residue_df.columns:
        raise KeyError("residue_df 必须包含 x 和 y 两列。")

    if center_mode == "xy_mean":
        center_x = float(residue_df["x"].mean())
        center_y = float(residue_df["y"].mean())

    elif center_mode == "xy_median":
        center_x = float(residue_df["x"].median())
        center_y = float(residue_df["y"].median())

    else:
        raise ValueError(
            f"不支持的 center_mode: {center_mode}\n"
            "当前仅支持 'xy_mean' 和 'xy_median'。"
        )

    return center_x, center_y


def add_radial_distance(
    residue_df: pd.DataFrame,
    center_x: float,
    center_y: float,
) -> pd.DataFrame:
    """
    为残基特征表添加径向距离 radial_distance。

    参数
    ----
    residue_df : pandas.DataFrame
        残基级特征表。
        至少需要包含：
        - x
        - y

    center_x : float
        孔道中心轴在 x 方向的位置。

    center_y : float
        孔道中心轴在 y 方向的位置。

    返回
    ----
    df : pandas.DataFrame
        添加 radial_distance 后的新表。

    计算公式
    ----
    radial_distance = sqrt((x - center_x)^2 + (y - center_y)^2)

    物理含义
    ----
    radial_distance 越小，说明残基 Cα 越靠近孔道中心轴。
    对于中空孔道结构，孔道内壁残基通常会位于相对靠近中心轴的位置。
    """

    df = residue_df.copy()

    required_cols = {"x", "y"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(f"residue_df 缺少必要列: {missing_cols}")

    df["center_x"] = center_x
    df["center_y"] = center_y

    df["radial_distance"] = np.sqrt(
        (df["x"] - center_x) ** 2 + (df["y"] - center_y) ** 2
    )

    return df


def add_xy_angle(
    residue_df: pd.DataFrame,
    center_x: float,
    center_y: float,
) -> pd.DataFrame:
    """
    为残基特征表添加 x-y 平面角度 theta。

    参数
    ----
    residue_df : pandas.DataFrame
        残基级特征表。

    center_x, center_y : float
        中心轴在 x-y 平面上的位置。

    返回
    ----
    df : pandas.DataFrame
        添加 theta_rad 和 theta_deg 后的新表。

    说明
    ----
    theta 表示残基围绕孔道中心轴的角度。
    对七聚体结构来说，这个角度以后可用于分析 7 个亚基在圆周方向上的分布。
    当前阶段它不是必须特征，但保留下来有利于后续分析。
    """

    df = residue_df.copy()

    dx = df["x"] - center_x
    dy = df["y"] - center_y

    # arctan2 返回弧度，范围是 [-pi, pi]
    df["theta_rad"] = np.arctan2(dy, dx)

    # 转换为角度，范围约为 [-180, 180]
    df["theta_deg"] = np.degrees(df["theta_rad"])

    return df


def add_z_normalized(residue_df: pd.DataFrame) -> pd.DataFrame:
    """
    添加 z 方向归一化坐标。

    参数
    ----
    residue_df : pandas.DataFrame
        残基级特征表。
        至少包含 z 列。

    返回
    ----
    df : pandas.DataFrame
        添加 z_min, z_max, z_norm 后的表格。

    z_norm 含义
    ----
    z_norm = (z - z_min) / (z_max - z_min)

    取值范围大致为 [0, 1]。

    为什么要加 z_norm？
    因为不同 PDB 文件的绝对坐标位置可能不同。
    归一化后的 z_norm 更适合不同结构之间进行粗略比较。
    """

    df = residue_df.copy()

    if "z" not in df.columns:
        raise KeyError("residue_df 必须包含 z 列。")

    z_min = float(df["z"].min())
    z_max = float(df["z"].max())

    df["z_min_global"] = z_min
    df["z_max_global"] = z_max

    if z_max == z_min:
        # 极端异常情况：所有 z 坐标相同
        df["z_norm"] = 0.0
    else:
        df["z_norm"] = (df["z"] - z_min) / (z_max - z_min)

    return df


def add_basic_geometry_features(
    residue_df: pd.DataFrame,
    center_mode: str = "xy_mean",
) -> pd.DataFrame:
    """
    为残基表添加基础几何特征。

    参数
    ----
    residue_df : pandas.DataFrame
        残基级特征表。

    center_mode : str
        孔道中心估计方式。
        默认 "xy_mean"。

    返回
    ----
    df : pandas.DataFrame
        添加如下列：
        - center_x
        - center_y
        - radial_distance
        - theta_rad
        - theta_deg
        - z_min_global
        - z_max_global
        - z_norm

    该函数是一个组合函数，用于让 pore.py 中的代码更简洁。
    """

    center_x, center_y = estimate_xy_center(
        residue_df=residue_df,
        center_mode=center_mode,
    )

    df = add_radial_distance(
        residue_df=residue_df,
        center_x=center_x,
        center_y=center_y,
    )

    df = add_xy_angle(
        residue_df=df,
        center_x=center_x,
        center_y=center_y,
    )

    df = add_z_normalized(df)

    return df