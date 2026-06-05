"""
01_sweep_spharm_bandwidth.py
============================
Bandwidth (h) SENSITIVITY SWEEP for the SP-SPHARM (scar-patterning) pipeline.

This is a NEW, self-contained add-on for the paper's Supplementary Information.
It does NOT modify the main pipeline, the cached results, or the manuscript.
It RE-USES the existing KDE / spherical-harmonic functions unchanged
(`batch_spherical_kde`, `kde_vector_to_dh_grid`, `compute_spharm_features`);
it only sweeps the von Mises-Fisher bandwidth h and writes the resulting
SP-SPHARM power spectra to a separate folder.

What it does, for each h in H_GRID:
    directions_aligned_svd.csv
        -> batch_spherical_kde(bandwidth = h)            (72 x 36 sphere grid)
        -> kde_vector_to_dh_grid(dh_size = 64)           (64 x 128 Driscoll-Healy)
        -> compute_spharm_features(lmax = 20)            (pyshtools expand + power)
        -> spectra/SPHARM_direction_h{h}.csv             (ID, Typology, n_scars,
                                                          power_l0 .. power_l20)

Only the SCAR (SP-SPHARM) side depends on h. Morphology (M-SPHARM) is
independent of h and is NOT recomputed here.

SANITY CHECK: at h = 0.35 the recomputed power spectrum must reproduce the
cached production file `derived_data/SPHARM_direction.csv`. The maximum absolute
per-degree difference is reported in `sweep_manifest.csv` and printed; a large
value (> 1e-6) is flagged loudly so it can be investigated (e.g. a BLAS / library
mismatch — see analysis/scripts/environment.yml, which pins the numerical core).

HOW TO RUN (canonical environment):
    # inside the project's conda `spharm` env (same one the main pipeline uses):
    python analysis/bandwidth_sensitivity/01_sweep_spharm_bandwidth.py

    # or via the Docker image, mirroring _targets.R's PYTHONPATH:
    #   PYTHONPATH=analysis/scripts python analysis/bandwidth_sensitivity/01_sweep_spharm_bandwidth.py

Outputs (all NEW, under analysis/bandwidth_sensitivity/):
    spectra/SPHARM_direction_h0.20.csv ... SPHARM_direction_h0.50.csv
    sweep_manifest.csv      (one row per h: kappa, n_specimens, max|diff| vs cache)
    versions.txt            (exact library versions used)

Then run 02_bandwidth_sensitivity_stats.R to evaluate stability of the
downstream conclusions across h.
"""

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import os
import sys
import platform
from pathlib import Path

# ---------------------------------------------------------------------------
# PARAMETERS  (edit here; coarsen H_GRID if the run is expensive)
# ---------------------------------------------------------------------------
# h grid. MUST include H_REF (0.35) for the sanity check.
# kappa = 1 / h^2, so this h range corresponds to kappa in [4.0, 25.0].
H_GRID = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
H_REF  = 0.35                       # the value used in the main analysis

# Fixed pipeline settings — kept IDENTICAL to kde_to_spharm_main.py so that the
# only thing that changes across the sweep is h.
N_BEARING = 72                      # KDE sphere-grid azimuth divisions
N_PLUNGE  = 36                      # KDE sphere-grid elevation divisions
LMAX      = 20                      # max spherical-harmonic degree
DH_SIZE   = 64                      # Driscoll-Healy latitude points (longitude = 2*DH_SIZE)
ALIGN_SRC = "svd"                   # production alignment (matches the main pipeline)

ROUND_DECIMALS = 8                  # match run_spharm() rounding in kde_to_spharm_main.py


# ---------------------------------------------------------------------------
# Locate the project root (the folder containing _targets.R) and wire up paths
# so the existing SPHARM_modules can be imported exactly as the main pipeline
# imports them.
# ---------------------------------------------------------------------------
def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "_targets.R").exists():
            return p
    # Fallback: this script lives at <root>/analysis/bandwidth_sensitivity/
    return start.parents[2]


THIS_DIR  = Path(__file__).resolve().parent
PROJ_ROOT = find_project_root(THIS_DIR)

SCRIPTS_DIR  = PROJ_ROOT / "analysis" / "scripts"
DERIVED_DIR  = PROJ_ROOT / "analysis" / "data" / "derived_data"
OUT_DIR      = THIS_DIR
SPECTRA_DIR  = OUT_DIR / "spectra"
INPUT_CSV    = DERIVED_DIR / f"directions_aligned_{ALIGN_SRC}.csv"
CACHE_CSV    = DERIVED_DIR / "SPHARM_direction.csv"     # cached h=0.35 production file

# Make the existing modules importable (mirrors _targets.R PYTHONPATH).
sys.path.insert(0, str(SCRIPTS_DIR))

# ---------------------------------------------------------------------------
# Third-party + REUSED project functions (NOT reimplemented)
# ---------------------------------------------------------------------------
import numpy as np
import pandas as pd

from SPHARM_modules.spherical_kde import batch_spherical_kde
from SPHARM_modules.kde_to_spharm import (
    kde_vector_to_dh_grid,
    compute_spharm_features,
)

# Determinism: the KDE -> SH -> power path has no stochastic component, but we
# seed anyway so the sweep is bit-for-bit identical to a single main-pipeline run.
np.random.seed(42)

POWER_COLS = [f"power_l{l}" for l in range(LMAX + 1)]


# ---------------------------------------------------------------------------
# One bandwidth -> power-spectrum table
# ---------------------------------------------------------------------------
def spectra_for_bandwidth(df: pd.DataFrame, h: float) -> pd.DataFrame:
    """
    Recompute the SP-SPHARM normalised power spectrum for every specimen at
    bandwidth h, reusing the production KDE + SH functions.

    Returns a DataFrame: ID, Typology, n_scars, power_l0 .. power_l{LMAX}.
    """
    # --- KDE (von Mises-Fisher spherical KDE), identical call to run_kde() ---
    kde_result = batch_spherical_kde(
        df,
        bandwidth    = h,
        n_bearing    = N_BEARING,
        n_plunge     = N_PLUNGE,
        id_col       = "ID",
        ux_col       = "ux",
        uy_col       = "uy",
        uz_col       = "uz",
        typology_col = "Typology",
        verbose      = False,
    )

    # --- sphere grid -> bearing/plunge frame, identical to run_spharm() ---
    G = kde_result["G"]
    sphere_grid = pd.DataFrame(G, columns=["x", "y", "z"])
    sphere_grid["bearing"] = np.arctan2(G[:, 1], G[:, 0])
    sphere_grid["plunge"]  = np.arcsin(np.clip(G[:, 2], -1, 1))

    kde_matrix = kde_result["kde_matrix"]
    rows = []
    for i, specimen_id in enumerate(kde_result["ids"]):
        typology = kde_result["typologies"][i]
        n_scars  = kde_result["n_scars"][i]
        try:
            grid_2d = kde_vector_to_dh_grid(kde_matrix[i], sphere_grid, dh_size=DH_SIZE)
            feats   = compute_spharm_features(grid_2d, lmax=LMAX)
            row = {"ID": specimen_id, "Typology": typology, "n_scars": n_scars}
            for l, p in enumerate(feats["norm_power"]):
                row[f"power_l{l}"] = round(float(p), ROUND_DECIMALS)
            rows.append(row)
        except Exception as e:                       # match run_spharm() robustness
            print(f"    [warn] {specimen_id}: {e}")
            rows.append({"ID": specimen_id, "Typology": typology, "n_scars": n_scars})

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Sanity check vs the cached production spectrum (only meaningful at h = H_REF)
# ---------------------------------------------------------------------------
def max_abs_diff_vs_cache(df_h: pd.DataFrame) -> float:
    if not CACHE_CSV.exists():
        print(f"    [info] cache not found ({CACHE_CSV.name}); skipping sanity check")
        return float("nan")
    cache = pd.read_csv(CACHE_CSV)
    cols  = [c for c in POWER_COLS if c in cache.columns and c in df_h.columns]
    a = df_h.set_index("ID")[cols].sort_index()
    b = cache.set_index("ID")[cols].reindex(a.index)
    return float(np.nanmax(np.abs(a.values - b.values)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 70)
    print("SP-SPHARM bandwidth sensitivity sweep")
    print(f"  project root : {PROJ_ROOT}")
    print(f"  input        : {INPUT_CSV}")
    print(f"  h grid       : {H_GRID}")
    print(f"  grid/lmax/DH : {N_BEARING}x{N_PLUNGE} / lmax={LMAX} / DH={DH_SIZE}x{2*DH_SIZE}")
    print("=" * 70)

    if not INPUT_CSV.exists():
        raise FileNotFoundError(
            f"Input not found: {INPUT_CSV}\n"
            "Run align_svd.R (targets target `align_svd_csvs`) first."
        )

    SPECTRA_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_CSV)
    missing = {"ID", "ux", "uy", "uz"} - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    if "Typology" not in df.columns:
        df["Typology"] = "unknown"
    print(f"Loaded {df['ID'].nunique()} specimens, {len(df)} scar vectors.\n")

    manifest = []
    for h in H_GRID:
        kappa = 1.0 / h ** 2
        print(f"[h = {h:.2f}]  kappa = {kappa:.3f}")
        df_h = spectra_for_bandwidth(df, h)

        out_csv = SPECTRA_DIR / f"SPHARM_direction_h{h:.2f}.csv"
        df_h.to_csv(out_csv, index=False)

        diff = float("nan")
        flag = ""
        if abs(h - H_REF) < 1e-9:
            diff = max_abs_diff_vs_cache(df_h)
            if not np.isnan(diff):
                ok = diff <= 1e-6
                flag = "OK" if ok else "*** MISMATCH ***"
                print(f"    sanity check vs cache: max|diff| = {diff:.3e}  [{flag}]")
                if not ok:
                    print("    >>> Recomputed h=0.35 spectra DIFFER from the cached "
                          "production file. Check library versions (numpy/pyshtools/"
                          "BLAS) against analysis/scripts/environment.yml.")

        manifest.append({
            "h"                 : h,
            "kappa"             : round(kappa, 6),
            "n_specimens"       : int(df_h["ID"].nunique()),
            "output_csv"        : str(out_csv.relative_to(PROJ_ROOT)),
            "is_reference"      : abs(h - H_REF) < 1e-9,
            "max_abs_diff_cache": diff,
            "sanity_flag"       : flag,
        })
        print(f"    wrote {out_csv.relative_to(PROJ_ROOT)}\n")

    pd.DataFrame(manifest).to_csv(OUT_DIR / "sweep_manifest.csv", index=False)
    print(f"Wrote {(OUT_DIR / 'sweep_manifest.csv').relative_to(PROJ_ROOT)}")

    # Record exact versions used (reproducibility).
    try:
        import scipy
        scipy_v = scipy.__version__
    except Exception:
        scipy_v = "not importable"
    try:
        import pyshtools
        pysh_v = pyshtools.__version__
    except Exception:
        pysh_v = "not importable"

    versions = [
        "SP-SPHARM bandwidth sweep — library versions",
        "=" * 46,
        f"timestamp        : {pd.Timestamp.now().isoformat()}",
        f"python           : {platform.python_version()} ({sys.executable})",
        f"platform         : {platform.platform()}",
        f"numpy            : {np.__version__}",
        f"pandas           : {pd.__version__}",
        f"scipy            : {scipy_v}",
        f"pyshtools        : {pysh_v}",
        "",
        "Pinned reference environment (analysis/scripts/environment.yml):",
        "  numpy=1.26.4, scipy=1.12.0, pyshtools=4.13.1, OpenBLAS 0.3.30",
        "If the versions above differ, the h=0.35 sanity check in sweep_manifest.csv",
        "confirms whether the numerical results still reproduce the cached spectra.",
    ]
    (OUT_DIR / "versions.txt").write_text("\n".join(versions), encoding="utf-8")
    print(f"Wrote {(OUT_DIR / 'versions.txt').relative_to(PROJ_ROOT)}")
    print("\nDone. Next: Rscript analysis/bandwidth_sensitivity/02_bandwidth_sensitivity_stats.R")


if __name__ == "__main__":
    main()
