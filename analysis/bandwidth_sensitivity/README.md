# SP-SPHARM bandwidth (h) sensitivity analysis

Supplementary, self-contained add-on that tests whether the paper's SP-SPHARM
(scar-patterning) results are robust to the von Mises–Fisher spherical-KDE
bandwidth, fixed at **h = 0.35** (κ = 1/h² ≈ 8.16) in the main analysis.

**This folder adds new files only.** It does not modify or overwrite the main
pipeline, the cached `_targets` store, the `derived_data` cache, or the
manuscript. The main results are unchanged.

## What depends on h

Only the **scar / SP-SPHARM** side depends on h. The KDE smooths the scar-vector
field before the spherical-harmonic expansion, so changing h changes the power
spectra `power_l*` in `SPHARM_direction.csv`. **Morphology (M-SPHARM) is
independent of h** and is reused unchanged (`SPHARM_morphology.csv` is never
recomputed).

## Files

| File | Role |
|------|------|
| `01_sweep_spharm_bandwidth.py` | Recomputes SP-SPHARM power spectra over the h grid, **reusing the existing KDE/SH functions** (`batch_spherical_kde`, `kde_vector_to_dh_grid`, `compute_spharm_features`). Writes `spectra/SPHARM_direction_h{h}.csv`. |
| `02_bandwidth_sensitivity_stats.R` | For each h, re-runs the downstream analyses with the **same package functions the main pipeline uses** (`vegan::adonis2`/`mantel`, `ade4::coinertia`/`randtest`, `compositions::ilr`), replicating the data-prep from the main scripts. Writes the metrics CSV, figures, and the SI summary. |
| `SI_methods_text.md` | Ready-to-paste SI **Methods** paragraph (no numbers needed). |
| *(generated)* `spectra/` | Per-h power-spectrum CSVs. |
| *(generated)* `sweep_manifest.csv`, `versions.txt` | Sweep bookkeeping + exact library versions. |
| *(generated)* `bandwidth_sensitivity_metrics.csv` | Tidy table, one row per h. |
| *(generated)* `bandwidth_orderselection_by_degree.csv` | Per-degree cumulative power & CV for each h. |
| *(generated)* `figures/fig_S_bandwidth_*.png` | Order-selection overlays, metric-vs-h lines, IM heatmaps. |
| *(generated)* `SI_bandwidth_sensitivity_summary.md` | Auto-filled SI **Results** summary + per-h table + robustness verdict. |

## How to run

Run inside the project's canonical environment (the same one the main pipeline
uses): the conda **`spharm`** env for Python and **R 4.4 + renv** for R — e.g.
via the project Docker image, or any environment that satisfies
`analysis/scripts/environment.yml` and `renv.lock`.

```bash
# 1) Recompute SP-SPHARM spectra across the h grid (reuses the production functions)
python analysis/bandwidth_sensitivity/01_sweep_spharm_bandwidth.py

# 2) Evaluate stability of the downstream conclusions and build SI outputs
Rscript analysis/bandwidth_sensitivity/02_bandwidth_sensitivity_stats.R
```

In Docker, mirror `_targets.R`'s settings (`PYTHONPATH=analysis/scripts`,
conda env `/opt/conda/envs/spharm`). Step 1 needs the production input
`analysis/data/derived_data/directions_aligned_svd.csv` (built by the
`align_svd_csvs` target); step 2 needs the committed `SPHARM_morphology.csv`
and the per-h spectra from step 1.

## Parameters

The h grid is a parameter at the top of **both** scripts — keep them in sync.
Coarsen it (e.g. `0.25, 0.35, 0.45`) if the run is expensive:

```python
H_GRID = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]   # 01_sweep_spharm_bandwidth.py
```
```r
H_GRID <- c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)  # 02_bandwidth_sensitivity_stats.R
```

Everything else is held identical to the main pipeline: SVD alignment, 72×36 KDE
grid, 64×128 Driscoll–Healy grid, l_max = 20, AC-normalisation, seed 42.

## Built-in reproducibility check

The grid **must include h = 0.35**. At that value the recomputed spectra must
reproduce the committed production file:

- `01_*.py` compares the recomputed h = 0.35 spectrum against
  `derived_data/SPHARM_direction.csv` and records `max|diff|` in
  `sweep_manifest.csv` (flagged if > 1e-6, e.g. a BLAS/library mismatch — the
  numerical core is pinned in `analysis/scripts/environment.yml`).
- `02_*.R` re-derives the h = 0.35 metrics and checks them against the committed
  values from `OrderSelection_stats_direction_EXP.csv`, `EXP_L1_results.csv`,
  `L1_results.csv`, and `L3_permanova.csv`. Deterministic statistics (pseudo-F,
  R², Mantel r, RV) must match; permutation p-values carry ~±0.005 Monte-Carlo
  jitter and are excluded from the check.

If the h = 0.35 anchor matches, the rest of the sweep is trustworthy and the
robustness statements in the SI summary follow directly.

## Metrics evaluated per h

- **(a) Order selection** — cumulative power through l = 6 and the degree at
  which cross-specimen CV first exceeds 100% (EXP and SDG): is the l = 1–6
  truncation robust?
- **(b) Ideal cores** — standardised-Euclidean distance structure (`power_l1:l4`)
  and its correlation with the h = 0.35 distance matrix; persistence of key
  separations (discoid vs Levallois, biface vs unifacial discoid).
- **(c) Experimental cores** — scar PERMANOVA core-type R²/pseudo-F/p and the set
  of significant pairwise distinctions (the "resolution profile").
- **(d) Decoupling** — global Mantel r and RV (EXP and SDG) and the SDG
  scar~core-type PERMANOVA: does the no-covariation conclusion hold across h?

## Optional: wiring into `targets`

The scripts are deliberately standalone so they touch no existing pipeline file.
If you later want them tracked by `targets`, add **new** targets to `_targets.R`
(do not edit existing ones) — e.g. a `format = "file"` target that runs
`01_*.py` and one that sources `02_*.R` — following the Python-invocation helpers
already defined there.
