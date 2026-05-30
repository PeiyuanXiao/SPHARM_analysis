"""
rotate_svd_directions.py
========================
Apply a random per-specimen Z-axis rotation to the SVD-aligned direction
vectors, producing directions_aligned_svd_rotated.csv for the empirical
rotational-invariance check.

Rationale:
  SVD alignment fixes the Z axis (the normal of the best-fit plane), but the
  in-plane (XY) rotation angle is arbitrary per specimen. Applying an
  independent random Z rotation to each specimen simulates the real
  directional uncertainty that exists within a fixed coordinate frame.

  If the SPHARM power spectrum is unchanged under this perturbation, then
  fixing the SVD frame makes the analysis insensitive to within-frame error.

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
RANDOM_SEED  = 42   # fixed seed for reproducibility


# =============================================================================
# Z-axis rotation matrix
# =============================================================================

def rot_z(theta: float) -> np.ndarray:
    """
    3x3 rotation matrix for a rotation of theta radians about the Z axis.
    The Z component is unchanged; rotation happens in the XY plane.
    """
    c, s = np.cos(theta), np.sin(theta)
    return np.array([
        [ c, -s, 0],
        [ s,  c, 0],
        [ 0,  0, 1],
    ])


# =============================================================================
# Main
# =============================================================================

def main():
    # --- load SVD-aligned data ---
    if not os.path.exists(INPUT_CSV):
        raise FileNotFoundError(
            f"Not found: {INPUT_CSV}\n"
            f"Run align_svd.R in R first to generate it."
        )

    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded: {INPUT_CSV}")
    print(f"  specimens: {df['ID'].nunique()}, scars: {len(df)}\n")

    # --- one independent random Z-rotation angle per specimen ---
    rng     = np.random.default_rng(RANDOM_SEED)
    all_ids = df["ID"].unique()

    # one random angle per specimen, in [0, 2*pi)
    angles  = rng.uniform(0, 2 * np.pi, size=len(all_ids))
    id_to_angle = dict(zip(all_ids, angles))

    print("Random rotation angle per specimen (radians):")
    for id_i, angle in id_to_angle.items():
        print(f"  {id_i:<40} {angle:.4f} rad  ({np.degrees(angle):.1f}°)")

    # --- apply rotation per specimen ---
    rotated_parts = []

    for id_i in all_ids:
        df_i   = df[df["ID"] == id_i].copy()
        theta  = id_to_angle[id_i]
        R      = rot_z(theta)

        # rotate direction vectors
        uv           = df_i[["ux", "uy", "uz"]].values  # (n, 3)
        uv_rotated   = uv @ R.T                          # (n, 3)

        df_i["ux"]   = uv_rotated[:, 0]
        df_i["uy"]   = uv_rotated[:, 1]
        df_i["uz"]   = uv_rotated[:, 2]

        rotated_parts.append(df_i)

    df_rotated = pd.concat(rotated_parts, ignore_index=True)

    # --- check: unit-vector length should stay 1 ---
    norms = np.sqrt(
        df_rotated["ux"]**2 +
        df_rotated["uy"]**2 +
        df_rotated["uz"]**2
    )
    print(f"\nUnit-vector length check (all should be ~1.0):")
    print(f"  min={norms.min():.8f}, max={norms.max():.8f}, "
          f"mean={norms.mean():.8f}")

    # --- save ---
    df_rotated.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved: {OUTPUT_CSV}")
    print(f"  specimens: {df_rotated['ID'].nunique()}, "
          f"scars: {len(df_rotated)}")


if __name__ == "__main__":
    main()
