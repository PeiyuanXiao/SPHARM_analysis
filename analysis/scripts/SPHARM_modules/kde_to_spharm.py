"""
kde_to_spharm.py
================
Spherical-harmonic helpers used by kde_to_spharm_main.py.

Contents:
    kde_vector_to_dh_grid()     - interpolate a KDE vector onto a DH grid
    compute_spharm_features()   - DH grid -> SH expansion -> power spectrum
    compute_variance_analysis() - per-degree power variance across specimens
"""

# Standard library
import sys
from pathlib import Path

# Path setup
sys.path.insert(0, str(Path(__file__).parent.parent))

# Third-party
import numpy as np
import pandas as pd
import pyshtools as pysh

# Local modules
from SPHARM_modules.power_spectrum import compute_power_spectrum

# ============================================================
# Default parameters
# ============================================================
LMAX    = 20
DH_SIZE = 64


# ============================================================
# Step 1: KDE vector -> standard 2D DH grid
# ============================================================
def kde_vector_to_dh_grid(kde_vec: np.ndarray,
                           sphere_grid: pd.DataFrame,
                           dh_size: int = DH_SIZE) -> np.ndarray:
    """
    Interpolate a KDE probability vector onto a Driscoll-Healy grid.

    Parameters
    ----------
    kde_vec     : ndarray, shape (n_grid,), KDE values for one specimen
    sphere_grid : DataFrame with columns bearing, plunge (radians)
    dh_size     : number of latitude points; longitude has 2*dh_size

    Returns
    -------
    grid_2d : ndarray, shape (dh_size, 2*dh_size), normalised DH grid
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
# Step 2: DH grid -> SH expansion -> power spectrum
# ============================================================
def compute_spharm_features(grid_2d: np.ndarray,
                             lmax: int = LMAX) -> dict:
    """
    Expand a DH grid in spherical harmonics and return its power spectrum.

    Parameters
    ----------
    grid_2d : ndarray, shape (dh_size, 2*dh_size)
    lmax    : maximum spherical-harmonic degree

    Returns
    -------
    dict with keys:
        power_spectrum - raw power spectrum
        norm_power     - normalised power spectrum
        coeffs_flat    - flattened spherical-harmonic coefficients
    """
    sh_grid = pysh.SHGrid.from_array(grid_2d, grid='DH')
    clm     = sh_grid.expand(lmax_calc=lmax)

    feats = compute_power_spectrum(clm, lmax=lmax)

    return {
        "power_spectrum" : feats["raw_power"],
        "norm_power"     : feats["norm_power"],
        "coeffs_flat"    : clm.coeffs.flatten(),
    }


# ============================================================
# Step 3: per-degree power variance across specimens
# ============================================================
def compute_variance_analysis(df_out: pd.DataFrame,
                               lmax: int,
                               output_csv: str) -> pd.DataFrame:
    """
    Per-degree variance of normalised power across all specimens.

    Parameters
    ----------
    df_out     : batch-result DataFrame with columns power_l0..power_lN
    lmax       : maximum spherical-harmonic degree
    output_csv : path to save the variance results

    Returns
    -------
    df_var : DataFrame with columns [degree, variance]
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
