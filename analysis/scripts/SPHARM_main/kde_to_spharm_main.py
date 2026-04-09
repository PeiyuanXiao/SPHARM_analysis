"""
kde_to_spharm_main.py
=====================
完整流水线：R 导出的方向向量 CSV → KDE → SPHARM → 功率谱 CSV

对齐步骤已移至 R 端（align_svd.R / align_lin2024.R），
本脚本读取 R 导出的方向向量，完成 KDE + SPHARM 两步。

两种使用模式：

  【日常生产】python kde_to_spharm_main.py
      默认 --source svd，输出到 derived_data/SPHARM_direction.csv
      与 R 脚本（spharm_analysis.R）的读取路径完全兼容

  【旋转不变性验证】python kde_to_spharm_main.py --source all
      依次处理 raw / svd / lin2024 三份数据
      各自输出到 derived_data/validation/{source}/SPHARM_direction.csv

Input CSVs (由 R 导出，位于 derived_data/):
    directions_raw.csv              — 原始数据（未对齐）
    directions_aligned_svd.csv      — SVD 法对齐
    directions_aligned_lin2024.csv  — Lin 2024 法对齐
"""

# ── 标准库 ──────────────────────────────────────
import os
import sys
import argparse
from pathlib import Path

# ── 路径设置 ─────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent.parent))  # → SPHARM_modules/

# ── 第三方库 ─────────────────────────────────────
import numpy as np
import pandas as pd

# ── 本地模块 ─────────────────────────────────────
from SPHARM_modules.spherical_kde import batch_spherical_kde
from SPHARM_modules.kde_to_spharm import (
    kde_vector_to_dh_grid,
    compute_spharm_features,
    compute_variance_analysis,
)


# =============================================================================
# Config
# =============================================================================

DERIVED_DIR = "/project/analysis/data/derived_data"
BANDWIDTH   = 0.35
N_BEARING   = 72
N_PLUNGE    = 36
LMAX        = 20
DH_SIZE     = 64

# R 导出的方向向量 CSV 路径
SOURCE_CSV = {
    "raw"         : f"{DERIVED_DIR}/directions_raw.csv",
    "svd"         : f"{DERIVED_DIR}/directions_aligned_svd.csv",
    "lin2024"     : f"{DERIVED_DIR}/directions_aligned_lin2024.csv",
    "svd_rotated" : f"{DERIVED_DIR}/directions_aligned_svd_rotated.csv",
}

# 日常生产模式的输出路径（svd 与 R 脚本兼容）
PRODUCTION_OUT = {
    "raw"         : f"{DERIVED_DIR}/SPHARM_direction_raw.csv",
    "svd"         : f"{DERIVED_DIR}/SPHARM_direction.csv",
    "lin2024"     : f"{DERIVED_DIR}/SPHARM_direction_lin2024.csv",
    "svd_rotated" : f"{DERIVED_DIR}/SPHARM_direction_svd_rotated.csv",
}

# 验证模式（--source all）的输出路径
VALIDATION_OUT = {
    "raw"         : f"{DERIVED_DIR}/validation/raw/SPHARM_direction.csv",
    "svd"         : f"{DERIVED_DIR}/validation/svd/SPHARM_direction.csv",
    "lin2024"     : f"{DERIVED_DIR}/validation/lin2024/SPHARM_direction.csv",
    "svd_rotated" : f"{DERIVED_DIR}/validation/svd_rotated/SPHARM_direction.csv",
}

# KDE 中间文件输出路径（供调试或下游脚本使用）
KDE_NPY_OUT = {
    "raw"         : f"{DERIVED_DIR}/kde_matrix_raw.npy",
    "svd"         : f"{DERIVED_DIR}/kde_matrix.npy",
    "lin2024"     : f"{DERIVED_DIR}/kde_matrix_lin2024.npy",
    "svd_rotated" : f"{DERIVED_DIR}/kde_matrix_svd_rotated.npy",
}


# =============================================================================
# Step 1: 读取 R 导出的方向向量 CSV
# =============================================================================

def load_directions(source: str) -> pd.DataFrame:
    csv_path = SOURCE_CSV[source]
    if not os.path.exists(csv_path):
        hint = {
            "raw"         : "请先在 R 中运行 align_svd.R",
            "svd"         : "请先在 R 中运行 align_svd.R",
            "lin2024"     : "请先在 R 中运行 align_lin2024.R",
            "svd_rotated" : "请先运行：python rotate_svd_directions.py",
        }.get(source, "请先生成该文件")
        raise FileNotFoundError(f"找不到：{csv_path}\n{hint}")

    df = pd.read_csv(csv_path)

    missing = {"ID", "ux", "uy", "uz"} - set(df.columns)
    if missing:
        raise ValueError(f"CSV 缺少必要列：{missing}")

    if "Typology" not in df.columns:
        df["Typology"] = "unknown"

    print(f"  Loaded: {df['ID'].nunique()} specimens, "
          f"{len(df)} direction vectors")
    return df


# =============================================================================
# Step 2: KDE，并保存中间文件
# =============================================================================

def run_kde(df: pd.DataFrame, source: str) -> dict:
    """
    对方向向量做球面 KDE，保存中间 .npy 文件供调试使用。

    Parameters
    ----------
    df     : 方向向量 DataFrame
    source : 数据源标识，用于确定 .npy 保存路径

    Returns
    -------
    kde_result : batch_spherical_kde 的完整返回字典
    """
    print(f"\n[KDE] bandwidth={BANDWIDTH}, grid={N_BEARING}×{N_PLUNGE}")
    kde_result = batch_spherical_kde(
        df,
        bandwidth    = BANDWIDTH,
        n_bearing    = N_BEARING,
        n_plunge     = N_PLUNGE,
        id_col       = "ID",
        ux_col       = "ux",
        uy_col       = "uy",
        uz_col       = "uz",
        typology_col = "Typology",
    )

    # 保存中间文件
    npy_path  = KDE_NPY_OUT[source]
    grid_path = npy_path.replace("kde_matrix", "kde_grid")
    meta_path = npy_path.replace("kde_matrix", "kde_metadata").replace(".npy", ".csv")

    np.save(npy_path, kde_result["kde_matrix"])
    np.save(grid_path, kde_result["G"])
    pd.DataFrame({
        "ID"      : kde_result["ids"],
        "Typology": kde_result["typologies"],
    }).to_csv(meta_path, index=False)

    print(f"  KDE 中间文件已保存：\n"
          f"    {npy_path}\n"
          f"    {grid_path}\n"
          f"    {meta_path}")

    return kde_result


# =============================================================================
# Step 3: KDE → DH 网格 → SPHARM，批量处理
# =============================================================================

def run_spharm(kde_result: dict,
               out_path: str,
               lmax: int = LMAX,
               dh_size: int = DH_SIZE) -> pd.DataFrame:
    """
    对 KDE 结果逐标本做球谐展开，保存功率谱 CSV 并输出方差分析。

    Parameters
    ----------
    kde_result : run_kde() 的返回值
    out_path   : 功率谱 CSV 保存路径
    lmax       : 最大球谐阶数
    dh_size    : DH 网格纬度方向点数

    Returns
    -------
    df_out : 功率谱 DataFrame
    """
    kde_matrix  = kde_result["kde_matrix"]
    sphere_grid = pd.DataFrame(kde_result["G"], columns=["x", "y", "z"])
    sphere_grid["bearing"] = np.arctan2(kde_result["G"][:, 1],
                                        kde_result["G"][:, 0])
    sphere_grid["plunge"]  = np.arcsin(np.clip(kde_result["G"][:, 2], -1, 1))

    print(f"\n[SPHARM] lmax={lmax}, DH grid={dh_size}×{dh_size*2}")
    rows = []

    for i, specimen_id in enumerate(kde_result["ids"]):
        typology = kde_result["typologies"][i]
        print(f"  [{i+1:>3}/{len(kde_result['ids'])}] {specimen_id}", end="  ")

        try:
            grid_2d = kde_vector_to_dh_grid(kde_matrix[i],
                                            sphere_grid,
                                            dh_size=dh_size)
            feats   = compute_spharm_features(grid_2d, lmax=lmax)

            row = {
                "ID"               : specimen_id,
                "Typology"         : typology,
                "spectral_entropy" : round(feats["spectral_entropy"], 6),
                "SHE"              : round(feats["she"], 6),
            }
            for l, p in enumerate(feats["norm_power"]):
                row[f"power_l{l}"] = round(float(p), 8)
            for j, c in enumerate(feats["coeffs_flat"]):
                row[f"coeff_{j}"] = round(float(np.real(c)), 8)

            rows.append(row)
            print(f"H={feats['spectral_entropy']:.4f}  ✓")

        except Exception as e:
            print(f"✗  {e}")
            rows.append({"ID": specimen_id, "Typology": typology})

    df_out = pd.DataFrame(rows)

    # 保存功率谱 CSV
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_out.to_csv(out_path, index=False)
    print(f"\n✓ Saved → {out_path}")
    print(f"  {len(df_out)} specimens, power_l0–power_l{lmax}")

    # 方差分析
    variance_csv = out_path.replace(".csv", "_variance_per_degree.csv")
    try:
        compute_variance_analysis(df_out, lmax, variance_csv)
    except Exception as e:
        print(f"Variance analysis failed (non-critical): {e}")

    return df_out


# =============================================================================
# 完整流水线
# =============================================================================

def run_pipeline(source: str,
                 validation: bool = False,
                 lmax: int = LMAX,
                 dh_size: int = DH_SIZE) -> pd.DataFrame:
    """
    对单个 source 运行完整流水线：
        方向向量 CSV → KDE → DH 网格 → SPHARM → 功率谱 CSV

    Parameters
    ----------
    source     : 'raw' | 'svd' | 'lin2024' | 'svd_rotated'
    validation : True  → 输出到 validation/ 子目录
                 False → 输出到生产路径（默认）
    """
    out_path = VALIDATION_OUT[source] if validation else PRODUCTION_OUT[source]

    print(f"\n{'='*60}")
    print(f"  Source : {source}")
    print(f"  Mode   : {'validation' if validation else 'production'}")
    print(f"  Output : {out_path}")
    print(f"{'='*60}")

    df         = load_directions(source)
    kde_result = run_kde(df, source)
    df_out     = run_spharm(kde_result, out_path, lmax=lmax, dh_size=dh_size)

    return df_out


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="KDE → SPHARM pipeline，读取 R 导出的方向向量 CSV"
    )
    parser.add_argument(
        "--source",
        choices=["raw", "svd", "lin2024", "svd_rotated", "all"],
        default="svd",
        help=(
            "raw         — 原始数据（未对齐）\n"
            "svd         — R SVD 法对齐（默认，日常生产使用）\n"
            "lin2024     — R Lin 2024 法对齐\n"
            "svd_rotated — SVD 对齐 + 随机 Z 轴旋转（实证验证用）\n"
            "all         — 依次运行四种，用于完整旋转不变性验证"
        )
    )
    args = parser.parse_args()

    validation_mode = (args.source == "all")
    sources = ["raw", "svd", "lin2024", "svd_rotated"] \
              if validation_mode else [args.source]

    for src in sources:
        run_pipeline(src, validation=validation_mode)

    print(f"\n{'='*60}")
    if validation_mode:
        print("验证模式完成。四组结果：")
        for src in sources:
            print(f"  {VALIDATION_OUT[src]}")
        print("\nNext: 运行 validate_rotation_all.R 比较各组功率谱")
    else:
        print(f"完成。输出：{PRODUCTION_OUT[args.source]}")
    print(f"{'='*60}")
