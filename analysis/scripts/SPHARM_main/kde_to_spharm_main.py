"""
kde_to_spharm_main.py
=====================
Full pipeline: R-exported direction-vector CSV -> KDE -> SPHARM -> power-spectrum CSV

Alignment is done in R (align_svd.R / align_lin2024.R);
this script reads the exported direction vectors and runs the KDE + SPHARM steps.

Two usage modes:

  [production] python kde_to_spharm_main.py
      defaults to --source svd, writing derived_data/SPHARM_direction.csv
      matches the path read by spharm_analysis.R

  [rotational-invariance validation] python kde_to_spharm_main.py --source all
      processes raw / svd / lin2024 in turn
      writing derived_data/validation/{source}/SPHARM_direction.csv each

Input CSVs (exported by R, in derived_data/):
    directions_raw.csv              - raw (unaligned)
    directions_aligned_svd.csv      - SVD alignment
    directions_aligned_lin2024.csv  - Lin 2024 alignment
"""

# Standard library
import os
import sys
import argparse
from pathlib import Path

# Path setup
sys.path.insert(0, str(Path(__file__).parent.parent))  # → SPHARM_modules/

# Third-party
import numpy as np
import pandas as pd

# Local modules
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

# R-exported direction-vector CSV paths
SOURCE_CSV = {
    "raw"         : f"{DERIVED_DIR}/directions_raw.csv",
    "svd"         : f"{DERIVED_DIR}/directions_aligned_svd.csv",
    "lin2024"     : f"{DERIVED_DIR}/directions_aligned_lin2024.csv",
    "svd_rotated" : f"{DERIVED_DIR}/directions_aligned_svd_rotated.csv",
}

# Production output paths (svd matches the R scripts)
PRODUCTION_OUT = {
    "raw"         : f"{DERIVED_DIR}/SPHARM_direction_raw.csv",
    "svd"         : f"{DERIVED_DIR}/SPHARM_direction.csv",
    "lin2024"     : f"{DERIVED_DIR}/SPHARM_direction_lin2024.csv",
    "svd_rotated" : f"{DERIVED_DIR}/SPHARM_direction_svd_rotated.csv",
}

# Validation-mode output paths (--source all)
VALIDATION_OUT = {
    "raw"         : f"{DERIVED_DIR}/validation/raw/SPHARM_direction.csv",
    "svd"         : f"{DERIVED_DIR}/validation/svd/SPHARM_direction.csv",
    "lin2024"     : f"{DERIVED_DIR}/validation/lin2024/SPHARM_direction.csv",
    "svd_rotated" : f"{DERIVED_DIR}/validation/svd_rotated/SPHARM_direction.csv",
}

# KDE intermediate-file paths (for debugging or downstream use)
KDE_NPY_OUT = {
    "raw"         : f"{DERIVED_DIR}/kde_matrix_raw.npy",
    "svd"         : f"{DERIVED_DIR}/kde_matrix.npy",
    "lin2024"     : f"{DERIVED_DIR}/kde_matrix_lin2024.npy",
    "svd_rotated" : f"{DERIVED_DIR}/kde_matrix_svd_rotated.npy",
}


# =============================================================================
# Step 1: load the R-exported direction-vector CSV
# =============================================================================

def load_directions(source: str) -> pd.DataFrame:
    csv_path = SOURCE_CSV[source]
    if not os.path.exists(csv_path):
        hint = {
            "raw"         : "Run align_svd.R in R first",
            "svd"         : "Run align_svd.R in R first",
            "lin2024"     : "Run align_lin2024.R in R first",
            "svd_rotated" : "Run python rotate_svd_directions.py first",
        }.get(source, "Generate this file first")
        raise FileNotFoundError(f"Not found: {csv_path}\n{hint}")

    df = pd.read_csv(csv_path)

    missing = {"ID", "ux", "uy", "uz"} - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    if "Typology" not in df.columns:
        df["Typology"] = "unknown"

    print(f"  Loaded: {df['ID'].nunique()} specimens, "
          f"{len(df)} direction vectors")
    return df


# =============================================================================
# Step 2: KDE, saving intermediate files
# =============================================================================

def run_kde(df: pd.DataFrame, source: str) -> dict:
    """
    Spherical KDE of the direction vectors; saves intermediate .npy files.

    Parameters
    ----------
    df     : direction-vector DataFrame
    source : data-source key, determines the .npy save path

    Returns
    -------
    kde_result : full dict returned by batch_spherical_kde
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

    # Save intermediate files
    npy_path  = KDE_NPY_OUT[source]
    grid_path = npy_path.replace("kde_matrix", "kde_grid")
    meta_path = npy_path.replace("kde_matrix", "kde_metadata").replace(".npy", ".csv")

    np.save(npy_path, kde_result["kde_matrix"])
    np.save(grid_path, kde_result["G"])
    pd.DataFrame({
        "ID"      : kde_result["ids"],
        "Typology": kde_result["typologies"],
    }).to_csv(meta_path, index=False)

    print(f"  KDE intermediate files saved:\n"
          f"    {npy_path}\n"
          f"    {grid_path}\n"
          f"    {meta_path}")

    return kde_result


# =============================================================================
# Step 3: KDE -> DH grid -> SPHARM, batch processing
# =============================================================================

def run_spharm(kde_result: dict,
               out_path: str,
               lmax: int = LMAX,
               dh_size: int = DH_SIZE) -> pd.DataFrame:
    """
    Per-specimen SH expansion of the KDE results; saves the power-spectrum CSV and variance analysis.

    Parameters
    ----------
    kde_result : return value of run_kde()
    out_path   : path to save the power-spectrum CSV
    lmax       : maximum spherical-harmonic degree
    dh_size    : number of DH-grid latitude points

    Returns
    -------
    df_out : power-spectrum DataFrame
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
                "ID"       : specimen_id,
                "Typology" : typology,
            }
            for l, p in enumerate(feats["norm_power"]):
                row[f"power_l{l}"] = round(float(p), 8)
            for j, c in enumerate(feats["coeffs_flat"]):
                row[f"coeff_{j}"] = round(float(np.real(c)), 8)

            rows.append(row)
            print("OK")

        except Exception as e:
            print(f"✗  {e}")
            rows.append({"ID": specimen_id, "Typology": typology})

    df_out = pd.DataFrame(rows)

    # Save the power-spectrum CSV
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_out.to_csv(out_path, index=False)
    print(f"\n✓ Saved → {out_path}")
    print(f"  {len(df_out)} specimens, power_l0–power_l{lmax}")

    # Variance analysis
    variance_csv = out_path.replace(".csv", "_variance_per_degree.csv")
    try:
        compute_variance_analysis(df_out, lmax, variance_csv)
    except Exception as e:
        print(f"Variance analysis failed (non-critical): {e}")

    return df_out


# =============================================================================
# Full pipeline
# =============================================================================

def run_pipeline(source: str,
                 validation: bool = False,
                 lmax: int = LMAX,
                 dh_size: int = DH_SIZE) -> pd.DataFrame:
    """
    Run the full pipeline for one source:
        direction-vector CSV -> KDE -> DH grid -> SPHARM -> power-spectrum CSV

    Parameters
    ----------
    source     : 'raw' | 'svd' | 'lin2024' | 'svd_rotated'
    validation : True  -> write to the validation/ subdirectory
                 False -> write to the production path (default)
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
        description="KDE -> SPHARM pipeline; reads the R-exported direction-vector CSV"
    )
    parser.add_argument(
        "--source",
        choices=["raw", "svd", "lin2024", "svd_rotated", "all"],
        default="svd",
        help=(
            "raw         - raw (unaligned)\n"
            "svd         - R SVD alignment (default, production)\n"
            "lin2024     - R Lin 2024 alignment\n"
            "svd_rotated - SVD alignment + random Z rotation (empirical validation)\n"
            "all         - run all four in turn, for full rotational-invariance validation"
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
        print("Validation mode complete. Four result sets:")
        for src in sources:
            print(f"  {VALIDATION_OUT[src]}")
        print("\nNext: run validate_rotation_all.R to compare the power spectra")
    else:
        print(f"Done. Output: {PRODUCTION_OUT[args.source]}")
    print(f"{'='*60}")
