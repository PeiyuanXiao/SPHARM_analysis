"""
========================
对 SVD 对齐后的方向向量施加随机 Z 轴旋转，
生成 directions_aligned_svd_rotated.csv 供旋转不变性实证验证使用。

设计逻辑：
  SVD 对齐已固定 Z 轴（法平面法线），但 XY 平面内的旋转角度
  对每个标本是任意的。对每个标本独立施加一个随机 Z 轴旋转角，
  模拟"固定坐标系内真实存在的方向不确定性"。

  如果 SPHARM 功率谱在此扰动下保持一致，
  说明固定 SVD 坐标系后，坐标系内的随机误差不影响后续分析。

Input:
  analysis/data/derived_data/directions_aligned_svd.csv

Output:
  analysis/data/derived_data/directions_aligned_svd_rotated.csv
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# =============================================================================
# Config
# =============================================================================

DERIVED_DIR  = "/project/analysis/data/derived_data"
INPUT_CSV    = f"{DERIVED_DIR}/directions_aligned_svd.csv"
OUTPUT_CSV   = f"{DERIVED_DIR}/directions_aligned_svd_rotated.csv"
RANDOM_SEED  = 42   # 固定种子保证可重复


# =============================================================================
# Z 轴旋转矩阵
# =============================================================================

def rot_z(theta: float) -> np.ndarray:
    """
    构造绕 Z 轴旋转 theta 弧度的 3×3 旋转矩阵。
    Z 轴分量不变，XY 平面内旋转。
    """
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [ c, -s, 0],
        [ s,  c, 0],
        [ 0,  0, 1],
    ])


# =============================================================================
# 主流程
# =============================================================================

def main():
    # --- 读取 SVD 对齐数据 ---
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(
            f"找不到：{INPUT_CSV}\n"
            f"请先在 R 中运行 align_svd.R 生成该文件。"
        )

    df = pd.read_csv(INPUT_CSV)
    print(f"读取：{INPUT_CSV}")
    print(f"  标本数：{df['ID'].nunique()}，刮痕总数：{len(df)}\n")

    # --- 对每个标本独立生成一个随机 Z 轴旋转角 ---
    rng     = np.random.default_rng(RANDOM_SEED)
    all_ids = df["ID"].unique()

    # 每个标本一个随机角度，范围 [0, 2π)
    angles  = rng.uniform(0, 2 * np.pi, size=len(all_ids))
    id_to_angle = dict(zip(all_ids, angles))

    print("各标本随机旋转角度（弧度）：")
    for id_i, angle in id_to_angle.items():
        print(f"  {id_i:<40} {angle:.4f} rad  ({np.degrees(angle):.1f}°)")

    # --- 逐标本施加旋转 ---
    rotated_parts = []

    for id_i in all_ids:
        df_i   = df[df["ID"] == id_i].copy()
        theta  = id_to_angle[id_i]
        R      = rot_z(theta)

        # 旋转方向向量
        uv           = df_i[["ux", "uy", "uz"]].values  # (n, 3)
        uv_rotated   = uv @ R.T                          # (n, 3)

        df_i["ux"]   = uv_rotated[:, 0]
        df_i["uy"]   = uv_rotated[:, 1]
        df_i["uz"]   = uv_rotated[:, 2]

        rotated_parts.append(df_i)

    df_rotated = pd.concat(rotated_parts, ignore_index=True)

    # --- 验证：单位向量长度应保持为 1 ---
    norms = np.sqrt(
        df_rotated["ux"]**2 +
        df_rotated["uy"]**2 +
        df_rotated["uz"]**2
    )
    print(f"\n单位向量长度检查（应全部 ≈ 1.0）：")
    print(f"  min={norms.min():.8f}, max={norms.max():.8f}, "
          f"mean={norms.mean():.8f}")

    # --- 保存 ---
    df_rotated.to_csv(OUTPUT_CSV, index=False)
    print(f"\n已保存：{OUTPUT_CSV}")
    print(f"  标本数：{df_rotated['ID'].nunique()}，"
          f"刮痕总数：{len(df_rotated)}")


if __name__ == "__main__":
    main()