"""
spharm_features.py
==================
球谐展开工具函数库，供 SPHARM_main/kde_to_spharm_main.py 调用。

包含：
    kde_vector_to_dh_grid()     — KDE 向量插值到 DH 标准网格
    compute_spharm_features()   — DH 网格 → 球谐展开 → 功率谱 + 谱熵
    compute_variance_analysis() — 跨标本每阶功率方差分析
"""

# ── 标准库 ──────────────────────────────────────
import sys
from pathlib import Path

# ── 路径设置 ─────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent))

# ── 第三方库 ─────────────────────────────────────
import numpy as np
import pandas as pd
import pyshtools as pysh

# ── 本地模块 ─────────────────────────────────────
from SPHARM_modules.spectral_entropy import compute_spectral_entropy

# ============================================================
# 默认参数
# ============================================================
LMAX    = 20
DH_SIZE = 64


# ============================================================
# Step 1: KDE 向量 → DH 标准二维网格
# ============================================================
def kde_vector_to_dh_grid(kde_vec: np.ndarray,
                           sphere_grid: pd.DataFrame,
                           dh_size: int = DH_SIZE) -> np.ndarray:
    """
    将 KDE 概率向量插值到 Driscoll-Healy 标准网格。

    Parameters
    ----------
    kde_vec     : ndarray, shape (n_grid,)，单个标本的 KDE 值
    sphere_grid : DataFrame，列包含 bearing、plunge（弧度）
    dh_size     : DH 网格纬度方向点数，经度方向为 2×dh_size

    Returns
    -------
    grid_2d : ndarray, shape (dh_size, 2*dh_size)，归一化后的 DH 网格
    """
    plunge  = sphere_grid["plunge"].values
    bearing = sphere_grid["bearing"].values

    colat_src = np.pi / 2 - plunge
    lon_src   = bearing

    n_lat    = dh_size
    n_lon    = 2 * dh_size
    colat_dh = np.linspace(0, np.pi,   n_lat, endpoint=False)
    lon_dh   = np.linspace(0, 2*np.pi, n_lon, endpoint=False)
    TH, PH   = np.meshgrid(colat_dh, lon_dh, indexing='ij')

    tx = np.sin(TH) * np.cos(PH)
    ty = np.sin(TH) * np.sin(PH)
    tz = np.cos(TH)

    sx = np.sin(colat_src) * np.cos(lon_src)
    sy = np.sin(colat_src) * np.sin(lon_src)
    sz = np.cos(colat_src)

    dot     = np.clip(
        tx[:, :, None]*sx + ty[:, :, None]*sy + tz[:, :, None]*sz,
        -1, 1
    )
    weights = np.exp(50 * dot)
    grid_2d = (np.sum(weights * kde_vec, axis=2) /
               np.sum(weights, axis=2))

    grid_2d = np.clip(grid_2d, 0, None)

    sin_weights = np.sin(colat_dh)[:, None]
    area_sum    = (grid_2d * sin_weights).sum()
    grid_2d    /= area_sum if area_sum > 0 else 1.0

    return grid_2d


# ============================================================
# Step 2: DH 网格 → 球谐展开 → 功率谱 + 谱熵
# ============================================================
def compute_spharm_features(grid_2d: np.ndarray,
                             lmax: int = LMAX) -> dict:
    """
    对 DH 网格做球谐展开，返回功率谱和谱熵。

    Parameters
    ----------
    grid_2d : ndarray, shape (dh_size, 2*dh_size)
    lmax    : 最大球谐阶数

    Returns
    -------
    dict，包含：
        power_spectrum   — 原始功率谱
        norm_power       — 归一化功率谱
        spectral_entropy — 谱熵 H
        she              — SHE 值
        coeffs_flat      — 展平的球谐系数
    """
    sh_grid = pysh.SHGrid.from_array(grid_2d, grid='DH')
    clm     = sh_grid.expand(lmax_calc=lmax)

    feats = compute_spectral_entropy(clm, lmax=lmax)

    return {
        "power_spectrum"   : feats["raw_power"],
        "norm_power"       : feats["norm_power"],
        "spectral_entropy" : feats["spectral_entropy"],
        "she"              : feats["SHE"],
        "coeffs_flat"      : clm.coeffs.flatten(),
    }


# ============================================================
# Step 3: 跨标本每阶功率方差分析
# ============================================================
def compute_variance_analysis(df_out: pd.DataFrame,
                               lmax: int,
                               output_csv: str) -> pd.DataFrame:
    """
    计算所有标本各球谐阶归一化功率的跨标本方差。

    Parameters
    ----------
    df_out     : batch 处理结果 DataFrame，含 power_l0–power_lN 列
    lmax       : 最大球谐阶数
    output_csv : 方差分析结果的保存路径

    Returns
    -------
    df_var : DataFrame，列为 [degree, variance]
    """
    power_cols = [f"power_l{l}" for l in range(lmax + 1)]
    available  = [c for c in power_cols if c in df_out.columns]

    variances = df_out[available].var(axis=0).values
    degrees   = list(range(len(variances)))

    df_var = pd.DataFrame({
        "degree":   degrees,
        "variance": variances,
    })

    df_var.to_csv(output_csv, index=False)

    print("\n==== Variance Analysis (Direction SPHARM) ====")
    print(f"Samples:       {len(df_out)}")
    print(f"Max variance:  {variances.max():.4f} (degree {variances.argmax()})")
    print(f"Min variance:  {variances.min():.4f} (degree {variances.argmin()})")
    print(f"Mean variance: {variances.mean():.4f}")
    print("\nTop 5 degrees by variance:")
    top5 = df_var.nlargest(5, "variance")
    for rank, (_, row) in enumerate(top5.iterrows(), 1):
        print(f"  Rank {rank}: degree {int(row['degree'])} → {row['variance']:.4f}")
    print(f"\nSaved variance analysis: {output_csv}")

    return df_var
