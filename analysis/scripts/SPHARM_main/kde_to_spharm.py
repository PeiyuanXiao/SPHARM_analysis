"""
kde_to_spharm.py
================
Convert spherical KDE outputs (CSV) from R into the standard DH grid format used by pyshtools,
and perform spherical harmonic expansion to generate power spectra, spherical harmonic coefficients, and spectral entropy.

Prerequisites: The spherical KDE must first be computed in R and the following files exported:
    output_dir <- "/project/analysis/data/derived_data"
    write.csv(sphere_grid, file.path(output_dir, "sphere_grid.csv"), row.names = FALSE)
    kde_df <- as.data.frame(kde_matrix)
    kde_df$ID <- rownames(kde_matrix)
    kde_df$Typology <- typology_vec
    write.csv(kde_df, file.path(output_dir, "kde_matrix.csv"), row.names = FALSE)

Input:
    kde_matrix.csv  - KDE value vector for each core (one row per specimen)
    sphere_grid.csv - Spherical grid coordinates (bearing, plunge, x, y, z)

Output:
    spharm_direction.csv - Power spectrum, spherical harmonic coefficients, and spectral entropy for each core
"""
import numpy as np
import pandas as pd
import pyshtools as pysh
from scipy.stats import entropy as scipy_entropy
 
 
# ============================================================
# Parameter configuration
# ============================================================
DATA_DIR = "/project/analysis/data/derived_data"
LMAX     = 20   
DH_SIZE  = 64  
 
 
# ============================================================
# Step 1: Load SKDE data exported from R
# ============================================================
def load_r_kde(
    kde_csv:  str = f"{DATA_DIR}/kde_matrix.csv",
    grid_csv: str = f"{DATA_DIR}/sphere_grid.csv"
):
    """
    Load kde_matrix.csv and sphere_grid.csv exported from R.

    Returns
    -------
    kde_matrix  : ndarray, shape (n_cores, n_grid)
    sphere_grid : DataFrame, containing columns: bearing, plunge, x, y, z
    meta        : DataFrame, containing ID and Typology
    """
    kde_df      = pd.read_csv(kde_csv)
    sphere_grid = pd.read_csv(grid_csv)

    meta_cols  = ["ID", "Typology"]
    meta       = kde_df[meta_cols].copy()
    kde_matrix = kde_df.drop(columns=meta_cols).values.astype(np.float64)

    print(f"Loaded {kde_matrix.shape[0]} cores, each with {kde_matrix.shape[1]} grid points")
    print(f"Typology distribution:\n{meta['Typology'].value_counts().to_string()}\n")
 
    return kde_matrix, sphere_grid, meta
 
 
# ============================================================
# Step 2: SKDE vector → DH standard 2D grid (vectorized vMF interpolation)
# ============================================================
def kde_vector_to_dh_grid(kde_vec: np.ndarray,
                           sphere_grid: pd.DataFrame,
                           dh_size: int = DH_SIZE) -> np.ndarray:
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
# Step 3: DH grid → spherical harmonic expansion → power spectrum + spectral entropy
# ============================================================
def compute_spharm_features(grid_2d: np.ndarray,
                             lmax: int = LMAX) -> dict:
    """
    Perform spherical harmonic expansion on the DH grid and return the power spectrum and spectral entropy
    """
    sh_grid = pysh.SHGrid.from_array(grid_2d, grid='DH')
    clm     = sh_grid.expand(lmax_calc=lmax)
 
    raw_power = clm.spectrum()          # shape: (lmax+1,)
    she       = float(np.sum(raw_power))
 
    # Normalize after removing the l = 0 DC component
    power_no_dc    = raw_power.copy()
    power_no_dc[0] = 0.0
    total_ac       = power_no_dc.sum()
    norm_power     = power_no_dc / total_ac if total_ac > 0 else power_no_dc
 
    # spectral entropy（l=1 ~ lmax）
    p = norm_power[1:]
    p = p[p > 0]
    spectral_entropy = float(scipy_entropy(p))
 
    return {
        "power_spectrum"   : raw_power,
        "norm_power"       : norm_power,
        "spectral_entropy" : spectral_entropy,
        "she"              : she,
        "coeffs_flat"      : clm.coeffs.flatten(),
    }
 
 
# ============================================================
# Step 4: Batch processing
# ============================================================
def batch_kde_to_spharm(
    kde_csv:    str = f"{DATA_DIR}/kde_matrix.csv",
    grid_csv:   str = f"{DATA_DIR}/sphere_grid.csv",
    output_csv: str = f"{DATA_DIR}/spharm_direction.csv",
    lmax:       int = LMAX,
    dh_size:    int = DH_SIZE
):
    """Batch process all cores and output a CSV containing power spectra and spectral entropy."""
    kde_matrix, sphere_grid, meta = load_r_kde(kde_csv, grid_csv)
    n_cores = kde_matrix.shape[0]
    rows    = []
 
    print(f"Starting spherical harmonic expansion (lmax={lmax}, DH grid={dh_size}×{dh_size*2})\n")
 
    for i in range(n_cores):
        specimen_id = meta["ID"].iloc[i]
        typology    = meta["Typology"].iloc[i]
        print(f"  [{i+1}/{n_cores}] {specimen_id}（{typology}）", end="  ")
 
        try:
            grid_2d = kde_vector_to_dh_grid(kde_matrix[i], sphere_grid, dh_size=dh_size)
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
            print(f"H={feats['spectral_entropy']:.4f}  SHE={feats['she']:.6f}  ✓")
 
        except Exception as e:
            print(f"✗ Error：{e}")
            rows.append({"ID": specimen_id, "Typology": typology})
 
    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_csv, index=False)
    print(f"\nCompleted! Results saved to: {output_csv}")
    print(f"Total {len(df_out)} cores processed; power spectrum columns: power_l0–power_l{lmax}, spectral entropy column: spectral_entropy")
 
    return df_out
 
 
# ============================================================
# Entry
# ============================================================
if __name__ == "__main__":
 
    df_results = batch_kde_to_spharm()
 
    print("\nMean spectral entropy by typology:")
    print(df_results.groupby("Typology")["spectral_entropy"].mean().round(4))
