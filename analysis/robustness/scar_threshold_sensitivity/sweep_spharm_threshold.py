"""sweep_spharm_threshold.py"""

import platform
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# PARAMETERS  (keep in sync with scar_attrition.py and 02_*.R)
# ---------------------------------------------------------------------------
THRESHOLDS = [0.0, 5.0, 10.0]
T_REF      = 0.0

BANDWIDTH = 0.35
N_BEARING = 72
N_PLUNGE  = 36
LMAX      = 20
DH_SIZE   = 64
ALIGN_SRC = "svd"

ROUND_DECIMALS = 8


# ---------------------------------------------------------------------------
# Locate project root (folder containing _targets.R) and wire up module imports
# exactly as the main pipeline does.
# ---------------------------------------------------------------------------
def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "_targets.R").exists():
            return p
    return start.parents[2]


THIS_DIR  = Path(__file__).resolve().parent
PROJ_ROOT = find_project_root(THIS_DIR)

SCRIPTS_DIR = PROJ_ROOT / "analysis" / "scripts"
RAW_DIR     = PROJ_ROOT / "analysis" / "data" / "raw_data"
DERIVED_DIR = PROJ_ROOT / "analysis" / "data" / "derived_data"
OUT_DIR     = THIS_DIR
SPECTRA_DIR = OUT_DIR / "spectra"

INPUT_CSV = DERIVED_DIR / f"directions_aligned_{ALIGN_SRC}.csv"
RAW_XLSX  = RAW_DIR / "Scar_orientation_data.xlsx"
CACHE_CSV = DERIVED_DIR / "SPHARM_direction.csv"

sys.path.insert(0, str(SCRIPTS_DIR))

import numpy as np
import pandas as pd

from SPHARM_modules.spherical_kde import batch_spherical_kde
from SPHARM_modules.kde_to_spharm import (
    kde_vector_to_dh_grid,
    compute_spharm_features,
)

np.random.seed(42)

POWER_COLS = [f"power_l{l}" for l in range(LMAX + 1)]


# ---------------------------------------------------------------------------
# Attach scar length (mm) to the production aligned vectors
# ---------------------------------------------------------------------------
def load_aligned_with_length() -> pd.DataFrame:
    """
    Return the production aligned direction vectors (ID, Typology, ux, uy, uz) with a
    `length_mm` column attached from the raw workbook, matched per specimen and
    per within-specimen scar order.
    """
    if not INPUT_CSV.exists():
        raise FileNotFoundError(
            f"Input not found: {INPUT_CSV}\n"
            "Run align_svd.R (targets target `align_svd_csvs`) first.")
    if not RAW_XLSX.exists():
        raise FileNotFoundError(f"Raw scar data not found: {RAW_XLSX}")

    aligned = pd.read_csv(INPUT_CSV)
    aligned["ID"] = aligned["ID"].astype(str).str.strip()
    missing = {"ID", "ux", "uy", "uz"} - set(aligned.columns)
    if missing:
        raise ValueError(f"Aligned CSV missing required columns: {missing}")
    if "Typology" not in aligned.columns:
        aligned["Typology"] = "unknown"

    xl = pd.ExcelFile(RAW_XLSX)
    frames = [xl.parse(xl.sheet_names[i]) for i in (0, 1, 2)]
    raw = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    raw["ID"] = raw["ID"].astype(str).str.strip()
    raw["length_mm"] = np.sqrt(
        (raw["End_X"] - raw["Start_X"]) ** 2
        + (raw["End_Y"] - raw["Start_Y"]) ** 2
        + (raw["End_Z"] - raw["Start_Z"]) ** 2)

    ca = aligned.groupby("ID").size()
    cr = raw.groupby("ID").size()
    bad = sorted(i for i in set(ca.index) | set(cr.index)
                 if ca.get(i, 0) != cr.get(i, 0))
    if bad:
        raise ValueError(
            "Per-specimen scar-count mismatch between aligned vectors and raw "
            f"workbook for: {bad[:10]}{'...' if len(bad) > 10 else ''}. "
            "Cannot re-attach scar lengths by position.")

    aligned = aligned.copy()
    aligned["seq"] = aligned.groupby("ID").cumcount()
    raw = raw.copy()
    raw["seq"] = raw.groupby("ID").cumcount()
    merged = aligned.merge(raw[["ID", "seq", "length_mm"]], on=["ID", "seq"], how="left")
    if merged["length_mm"].isna().any():
        n = int(merged["length_mm"].isna().sum())
        raise ValueError(f"{n} aligned scars could not be matched to a raw length.")
    return merged


def assemblage_of(idstr: str) -> str:
    s = str(idstr).strip()
    return "IM" if s.startswith("IM_") else ("SDG" if s.startswith("SDG")
            else ("EXP" if s.startswith("EXP") else "OTHER"))


def filter_for_threshold(df: pd.DataFrame, T: float) -> pd.DataFrame:
    """Keep a scar if it belongs to a synthetic IM core (held fixed) or length > T."""
    is_im = df["ID"].str.startswith("IM_")
    keep  = is_im | (df["length_mm"] > T)
    return df.loc[keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# One threshold -> power-spectrum table  (identical KDE/SH calls to run_spharm)
# ---------------------------------------------------------------------------
def spectra_for_threshold(df: pd.DataFrame) -> pd.DataFrame:
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
        verbose      = False,
    )

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
        except Exception as e:
            print(f"    [warn] {specimen_id}: {e}")
            rows.append({"ID": specimen_id, "Typology": typology, "n_scars": n_scars})

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Sanity check vs the cached production spectrum (only meaningful at T = T_REF)
# ---------------------------------------------------------------------------
def max_abs_diff_vs_cache(df_t: pd.DataFrame) -> float:
    if not CACHE_CSV.exists():
        print(f"    [info] cache not found ({CACHE_CSV.name}); skipping sanity check")
        return float("nan")
    cache = pd.read_csv(CACHE_CSV)
    cache["ID"] = cache["ID"].astype(str).str.strip()
    cols  = [c for c in POWER_COLS if c in cache.columns and c in df_t.columns]
    a = df_t.assign(ID=df_t["ID"].astype(str).str.strip()).set_index("ID")[cols].sort_index()
    b = cache.set_index("ID")[cols].reindex(a.index)
    return float(np.nanmax(np.abs(a.values - b.values)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 70)
    print("SP-SPHARM scar minimum-size-threshold sensitivity sweep")
    print(f"  project root : {PROJ_ROOT}")
    print(f"  aligned in   : {INPUT_CSV.relative_to(PROJ_ROOT)}")
    print(f"  lengths from : {RAW_XLSX.relative_to(PROJ_ROOT)}")
    print(f"  thresholds   : {THRESHOLDS} mm  (keep length > T; T=0 = all scars)")
    print(f"  h/grid/lmax  : {BANDWIDTH} / {N_BEARING}x{N_PLUNGE} / lmax={LMAX} / DH={DH_SIZE}x{2*DH_SIZE}")
    print("=" * 70)

    SPECTRA_DIR.mkdir(parents=True, exist_ok=True)

    df = load_aligned_with_length()
    df["group"] = df["ID"].map(assemblage_of)
    print(f"Loaded {df['ID'].nunique()} specimens, {len(df)} scars "
          f"(EXP={ (df.group=='EXP').sum() }, SDG={ (df.group=='SDG').sum() }, "
          f"IM={ (df.group=='IM').sum() } scars).\n")

    manifest = []
    for T in THRESHOLDS:
        df_t = filter_for_threshold(df, T)
        g = df_t["ID"].map(assemblage_of)
        n_exp, n_sdg, n_im = int((g == "EXP").sum()), int((g == "SDG").sum()), int((g == "IM").sum())
        print(f"[T = {T:>4.1f} mm]  scars kept: EXP={n_exp}, SDG={n_sdg}, IM={n_im} (held), "
              f"total={len(df_t)}; specimens={df_t['ID'].nunique()}")

        spec = spectra_for_threshold(df_t)
        out_csv = SPECTRA_DIR / f"SPHARM_direction_t{T:04.1f}.csv"
        spec.to_csv(out_csv, index=False)

        diff, flag = float("nan"), ""
        if abs(T - T_REF) < 1e-9:
            diff = max_abs_diff_vs_cache(spec)
            if not np.isnan(diff):
                ok = diff <= 1e-6
                flag = "OK" if ok else "*** MISMATCH ***"
                print(f"    sanity check vs cache (T=0): max|diff| = {diff:.3e}  [{flag}]")
                if not ok:
                    print("    >>> Recomputed T=0 spectra DIFFER from the cached production "
                          "file. Check library versions (numpy/pyshtools/BLAS) against "
                          "analysis/scripts/environment.yml.")

        manifest.append({
            "threshold_mm"      : T,
            "is_reference"      : abs(T - T_REF) < 1e-9,
            "n_specimens"       : int(spec["ID"].nunique()),
            "scars_kept_total"  : int(len(df_t)),
            "scars_kept_EXP"    : n_exp,
            "scars_kept_SDG"    : n_sdg,
            "scars_kept_IM"     : n_im,
            "output_csv"        : str(out_csv.relative_to(PROJ_ROOT)),
            "max_abs_diff_cache": diff,
            "sanity_flag"       : flag,
        })
        print(f"    wrote {out_csv.relative_to(PROJ_ROOT)}\n")

    pd.DataFrame(manifest).to_csv(OUT_DIR / "sweep_manifest.csv", index=False)
    print(f"Wrote {(OUT_DIR / 'sweep_manifest.csv').relative_to(PROJ_ROOT)}")

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

    print("\nDone. Next: Rscript analysis/scar_threshold_sensitivity/"
          "scar_threshold_sensitivity_stats.R")


if __name__ == "__main__":
    main()
