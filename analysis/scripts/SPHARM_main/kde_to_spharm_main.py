"""
kde_to_spharm_main.py
=====================
Pipeline:
    1. Load raw scar orientation data from Excel
    2. Align direction vectors (SVD normal → Z-axis, PCA main axis → X-axis)
    3. Run spherical KDE on aligned unit vectors
    4. Export results

Output files (all in /project/analysis/data/derived_data/):
    kde_matrix.npy      — (n_specimens, n_grid) density matrix
    kde_grid.npy        — (n_grid, 3) sphere grid unit vectors
    kde_metadata.csv    — specimen ID, typology, n_scars
    kde_results.csv     — long-format KDE densities (optional, for R)
"""

import os
import sys
import numpy as np
import pandas as pd

# Allow importing from SPHARM_modules
SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, SCRIPTS_DIR)

from SPHARM_modules.spherical_kde import batch_spherical_kde

# =============================================================================
# Config
# =============================================================================

INPUT_XLSX  = "/project/analysis/data/raw_data/Scar_orientation_data.xlsx"
OUTPUT_DIR  = "/project/analysis/data/derived_data"
BANDWIDTH   = 0.35
N_BEARING   = 72
N_PLUNGE    = 36


# =============================================================================
# Step 1: Load data
# =============================================================================

def load_data(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    print(f"Loaded {len(df)} rows from {os.path.basename(path)}")
    print(f"Columns: {list(df.columns)}\n")
    return df


# =============================================================================
# Step 2: Alignment (SVD normal → Z, PCA main axis → X)
# =============================================================================

def get_rot_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rodrigues rotation matrix that rotates unit vector a onto unit vector b."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    cos_theta = np.clip(np.dot(a, b), -1.0, 1.0)

    if cos_theta < -1 + 1e-10:
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        v = perp - np.dot(perp, a) * a
        v /= np.linalg.norm(v)
        return 2 * np.outer(v, v) - np.eye(3)

    if cos_theta > 1 - 1e-10:
        return np.eye(3)

    v = np.cross(a, b)
    vx = np.array([
        [ 0,    -v[2],  v[1]],
        [ v[2],  0,    -v[0]],
        [-v[1],  v[0],  0   ]
    ])
    return np.eye(3) + vx + vx @ vx * ((1 - cos_theta) / np.dot(v, v))


def align_group(df_group: pd.DataFrame) -> pd.DataFrame:
    """
    Align one specimen:
      Step 1 — SVD normal → Z-axis
      Step 2 — Translate plane centre to origin
      Step 3 — SVD main axis in XY plane → X-axis
    """
    df = df_group.copy()

    # --- Direction vectors ---
    dx  = df["End_X"].values - df["Start_X"].values
    dy  = df["End_Y"].values - df["Start_Y"].values
    dz  = df["End_Z"].values - df["Start_Z"].values
    length = np.sqrt(dx**2 + dy**2 + dz**2)
    valid  = length > 1e-10

    if valid.sum() < 3:
        print(f"  Warning: {df['ID'].iloc[0]} has fewer than 3 valid scars, skipping.")
        return None

    # Step 1: SVD normal → Z-axis
    U = np.column_stack([
        dx[valid] / length[valid],
        dy[valid] / length[valid],
        dz[valid] / length[valid],
    ])
    normal = np.linalg.svd(U)[2][-1]           # last right singular vector
    normal = normal / np.linalg.norm(normal)
    if normal[2] < 0:
        normal = -normal

    R1 = get_rot_matrix(normal, np.array([0, 0, 1]))

    S1 = df[["Start_X", "Start_Y", "Start_Z"]].values @ R1.T
    E1 = df[["End_X",   "End_Y",   "End_Z"  ]].values @ R1.T
    D1 = df[["Direct_X","Direct_Y","Direct_Z"]].values @ R1.T

    # Step 2: Translate centre to origin
    centre    = df[["Pos_X", "Pos_Y", "Pos_Z"]].iloc[0].values
    centre_r1 = R1 @ centre

    S1 -= centre_r1
    E1 -= centre_r1

    # Step 3: SVD main axis in XY → X-axis
    d2     = E1 - S1
    len2   = np.linalg.norm(d2, axis=1)
    valid2 = len2 > 1e-10

    xy_dirs  = d2[valid2, :2]
    main_dir = np.linalg.svd(xy_dirs)[2][0]     # first right singular vector
    mean_dir = xy_dirs.mean(axis=0)
    if np.dot(main_dir, mean_dir) < 0:
        main_dir = -main_dir

    theta = np.arctan2(main_dir[1], main_dir[0])
    c, s  = np.cos(-theta), np.sin(-theta)
    R2    = np.array([
        [ c, -s, 0],
        [ s,  c, 0],
        [ 0,  0, 1],
    ])

    S2 = S1 @ R2.T
    E2 = E1 @ R2.T
    D2 = D1 @ R2.T

    # Write aligned columns
    df[["s_x","s_y","s_z"]] = S2
    df[["e_x","e_y","e_z"]] = E2
    df[["d_x","d_y","d_z"]] = D2

    # Normalised direction unit vectors for KDE
    dv   = E2 - S2
    norm = np.linalg.norm(dv, axis=1, keepdims=True)
    norm = np.where(norm < 1e-10, 1.0, norm)
    uv   = dv / norm
    df["ux"] = uv[:, 0]
    df["uy"] = uv[:, 1]
    df["uz"] = uv[:, 2]

    return df


def align_all(raw: pd.DataFrame) -> pd.DataFrame:
    print("Aligning specimens...")
    aligned_parts = []
    for id_i, grp in raw.groupby("ID", sort=False):
        result = align_group(grp)
        if result is not None:
            aligned_parts.append(result)
    aligned = pd.concat(aligned_parts, ignore_index=True)
    print(f"  {aligned['ID'].nunique()} specimens aligned, "
          f"{len(aligned)} scars total.\n")
    return aligned


# =============================================================================
# Step 3 & 4: KDE + Export
# =============================================================================

def export_results(kde_result: dict, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    # Binary arrays (fast, used by Python downstream)
    np.save(os.path.join(output_dir, "kde_matrix.npy"), kde_result["kde_matrix"])
    np.save(os.path.join(output_dir, "kde_grid.npy"),   kde_result["G"])
    print(f"  Saved kde_matrix.npy  {kde_result['kde_matrix'].shape}")
    print(f"  Saved kde_grid.npy    {kde_result['G'].shape}")

    # Metadata CSV
    meta = pd.DataFrame({
        "ID":       kde_result["ids"],
        "Typology": kde_result["typologies"],
        "n_scars":  kde_result["n_scars"],
    })
    meta_path = os.path.join(output_dir, "kde_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"  Saved kde_metadata.csv  ({len(meta)} rows)")

    # Long-format CSV (optional, for R / manual inspection)
    rows = []
    for i, id_i in enumerate(kde_result["ids"]):
        for j, dens in enumerate(kde_result["kde_matrix"][i]):
            rows.append({"ID": id_i, "grid_point": j, "density": dens})
    long_df   = pd.DataFrame(rows)
    long_path = os.path.join(output_dir, "kde_results_long.csv")
    long_df.to_csv(long_path, index=False)
    print(f"  Saved kde_results_long.csv  ({len(long_df)} rows)")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Spherical KDE Pipeline")
    print("=" * 60)

    raw     = load_data(INPUT_XLSX)
    aligned = align_all(raw)

    kde_result = batch_spherical_kde(
        aligned,
        bandwidth = BANDWIDTH,
        n_bearing = N_BEARING,
        n_plunge  = N_PLUNGE,
    )

    print("Exporting results...")
    export_results(kde_result, OUTPUT_DIR)

    print("\n" + "=" * 60)
    print("Done!")
    print(f"  Specimens : {len(kde_result['ids'])}")
    print(f"  Grid size : {len(kde_result['G'])} points")
    print(f"  Output    : {OUTPUT_DIR}")
    print("=" * 60)