"""
sweep_spharm_threshold.py
============================
Scar minimum-SIZE-THRESHOLD sensitivity sweep for the SP-SPHARM (scar-patterning)
pipeline. NEW, self-contained add-on for the paper's Supplementary Information.

It does NOT modify the main pipeline, the cached results, or the manuscript, and it
RE-USES the existing KDE / spherical-harmonic functions unchanged
(`batch_spherical_kde`, `kde_vector_to_dh_grid`, `compute_spharm_features`). The
only thing that changes across the sweep is which scars are included in the KDE,
controlled by a minimum-size cutoff on scar length.

WHAT IT DOES, for each threshold T in THRESHOLDS:
    directions_aligned_svd.csv  (production aligned unit vectors)
        + scar length (mm) re-attached from Scar_orientation_data.xlsx
        -> drop scars with length <= T  (EXP & SDG only; IM held — see below)
        -> batch_spherical_kde(bandwidth = 0.35)        (72 x 36 sphere grid)
        -> kde_vector_to_dh_grid(dh_size = 64)          (64 x 128 Driscoll-Healy)
        -> compute_spharm_features(lmax = 20)           (pyshtools expand + power)
        -> spectra/SPHARM_direction_t{T}.csv            (ID, Typology, n_scars,
                                                         power_l0 .. power_l20)

DESIGN DECISIONS (documented; see README.md for the full pipeline map)
----------------------------------------------------------------------
* Scar size = the 3D Euclidean Start->End length (mm). This is exactly the `len`
  align_svd.R already computes, and it is rotation-invariant, so the size of a scar
  is unchanged by the SVD alignment. We attach it back to the *production* aligned
  vectors so the ALIGNMENT IS HELD FIXED at the published solution: the sweep
  isolates the effect of dropping small scars from the KDE, not a re-fitted
  alignment. (We do not re-run the SVD plane fit on the reduced scar sets.)

* The cutoff uses strict ">" (keep length > T), matching the manuscript's wording
  ("only scars larger than 5 mm"). With continuous lengths, >= vs > differs only
  for exact-boundary scars (none occur here).

* The production pipeline applies NO size filter (align_svd.R keeps len > 1e-10), so
  the committed spectra correspond to T = 0 (all recorded scars, realised min
  ~2 mm). T = 0 is therefore the REPRODUCIBILITY ANCHOR; 5 and 10 mm are the
  reviewer-facing perturbations. The manuscript's stated ">5 mm" is an acquisition
  guideline that is NOT enforced downstream — flagged here so it is not mistaken
  for the anchor.

* IDEAL (IM) cores are held at production values and are never re-filtered: their
  scar lengths are synthetic (IM_discoid / IM_discoid_unifacial are uniformly
  ~2.1 mm and would be erased by any >=5 mm cut). Their SP-SPHARM is unchanged by
  construction. The keep-rule is therefore: keep a scar if its ID is an IM_ core OR
  its length > T.

SANITY CHECK: at T = 0 the recomputed power spectrum must reproduce the cached
production file `derived_data/SPHARM_direction.csv`. The maximum absolute per-degree
difference is recorded in `sweep_manifest.csv` and printed; > 1e-6 is flagged (e.g.
a BLAS/library mismatch — the numerical core is pinned in
analysis/scripts/environment.yml).

HOW TO RUN (canonical environment — conda `spharm`, same as the main pipeline):
    python analysis/scar_threshold_sensitivity/sweep_spharm_threshold.py
    # or, mirroring _targets.R: PYTHONPATH=analysis/scripts python .../sweep_spharm_threshold.py

Then run scar_threshold_sensitivity_stats.R to evaluate stability of the
downstream conclusions across thresholds.
"""

import platform
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# PARAMETERS  (keep in sync with scar_attrition.py and 02_*.R)
# ---------------------------------------------------------------------------
# Minimum-size cutoffs in mm. MUST include T_REF (0.0) for the sanity check.
# Keep scars with length > T. T = 0 reproduces the production (all-scar) spectra.
THRESHOLDS = [0.0, 5.0, 10.0]
T_REF      = 0.0                    # production anchor (no size filter in the main pipeline)

# Fixed pipeline settings — kept IDENTICAL to kde_to_spharm_main.py so the only
# thing that changes across the sweep is the scar inclusion set.
BANDWIDTH = 0.35                    # vMF bandwidth used in the main analysis
N_BEARING = 72                      # KDE sphere-grid azimuth divisions
N_PLUNGE  = 36                      # KDE sphere-grid elevation divisions
LMAX      = 20                      # max spherical-harmonic degree
DH_SIZE   = 64                      # Driscoll-Healy latitude points (longitude = 2*DH_SIZE)
ALIGN_SRC = "svd"                   # production alignment (matches the main pipeline)

ROUND_DECIMALS = 8                  # match run_spharm() rounding in kde_to_spharm_main.py


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

INPUT_CSV = DERIVED_DIR / f"directions_aligned_{ALIGN_SRC}.csv"     # production aligned vectors
RAW_XLSX  = RAW_DIR / "Scar_orientation_data.xlsx"                  # source of scar lengths
CACHE_CSV = DERIVED_DIR / "SPHARM_direction.csv"                    # cached T=0 production file

sys.path.insert(0, str(SCRIPTS_DIR))     # make SPHARM_modules importable (mirrors _targets.R)

import numpy as np
import pandas as pd

from SPHARM_modules.spherical_kde import batch_spherical_kde
from SPHARM_modules.kde_to_spharm import (
    kde_vector_to_dh_grid,
    compute_spharm_features,
)

# Determinism: the KDE -> SH -> power path is deterministic, but seed anyway so the
# sweep is bit-for-bit identical to a single main-pipeline run.
np.random.seed(42)

POWER_COLS = [f"power_l{l}" for l in range(LMAX + 1)]


# ---------------------------------------------------------------------------
# Attach scar length (mm) to the production aligned vectors
# ---------------------------------------------------------------------------
def load_aligned_with_length() -> pd.DataFrame:
    """
    Return the production aligned direction vectors (ID, Typology, ux, uy, uz) with a
    `length_mm` column attached from the raw workbook, matched per specimen and
    per within-specimen scar order. Alignment is a pure rotation/translation, so the
    scar length is identical before and after alignment; we therefore re-attach the
    raw length to each production row by (ID, within-ID position).

    Robustness: IDs are stripped of stray whitespace (the workbook contains e.g.
    "IM_Multiplatform "), and we assert that every specimen has the same scar count
    in both files before merging.
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

    # Raw scars (sheets 1-3, in the same order align_svd.R binds them).
    xl = pd.ExcelFile(RAW_XLSX)
    frames = [xl.parse(xl.sheet_names[i]) for i in (0, 1, 2)]
    raw = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    raw["ID"] = raw["ID"].astype(str).str.strip()
    raw["length_mm"] = np.sqrt(
        (raw["End_X"] - raw["Start_X"]) ** 2
        + (raw["End_Y"] - raw["Start_Y"]) ** 2
        + (raw["End_Z"] - raw["Start_Z"]) ** 2)

    # Per-specimen count check (after stripping IDs).
    ca = aligned.groupby("ID").size()
    cr = raw.groupby("ID").size()
    bad = sorted(i for i in set(ca.index) | set(cr.index)
                 if ca.get(i, 0) != cr.get(i, 0))
    if bad:
        raise ValueError(
            "Per-specimen scar-count mismatch between aligned vectors and raw "
            f"workbook for: {bad[:10]}{'...' if len(bad) > 10 else ''}. "
            "Cannot re-attach scar lengths by position.")

    # Within-ID running index in each file's natural order (preserved by align_svd.R's
    # group_modify), then merge on (ID, seq).
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
        except Exception as e:                        # match run_spharm() robustness
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
        # per-group retained scar counts (cross-check against scar_attrition.py)
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

    print("\nDone. Next: Rscript analysis/scar_threshold_sensitivity/"
          "scar_threshold_sensitivity_stats.R")


if __name__ == "__main__":
    main()
