"""
00_mesh_inventory.py
====================
Mesh INVENTORY / decimation-context report for the mesh-preprocessing sensitivity
analysis (SI add-on). Self-contained; reads only committed data and writes new
files under analysis/mesh_decimation_sensitivity/. Does NOT touch the main pipeline.

WHY THIS IS A SEPARATE, LIGHTWEIGHT SCRIPT
------------------------------------------
The M-SPHARM re-decimation sweep (01_*.py) needs the project's meshing stack
(open3d + trimesh + pyshtools, the conda `spharm` env / Docker). But the *context*
for the sweep — how many faces each raw mesh has, and therefore which specimens a
given decimation target actually changes — depends only on reading the STL headers
and the committed original face counts, so it is computed here with the standard
library + pandas alone and runs on any machine.

KEY POINT THIS REPORT MAKES EXPLICIT
------------------------------------
`open3d.simplify_quadric_decimation(TARGET_FACES)` only ever REDUCES the face
count. A mesh whose raw face count is already below a given target is therefore
returned unchanged at that target, so its M-SPHARM spectrum is identical to
production there. This report counts, for each candidate target in FACE_TARGETS,
how many specimens are actually decimated vs left untouched — which is exactly what
makes a face-count setting "bite" downstream.

OUTPUTS (all NEW, under analysis/mesh_decimation_sensitivity/)
    mesh_inventory.csv          per specimen: raw faces, format, file size, group
    mesh_inventory_summary.csv  per group x target: # decimated vs unchanged
    figures/fig_S_mesh_inventory.png   (best-effort; needs matplotlib)

HOW TO RUN
    python analysis/mesh_decimation_sensitivity/00_mesh_inventory.py
"""

from __future__ import annotations

import platform
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# PARAMETERS  (keep FACE_TARGETS / production in sync with 01_*.py and 02_*.R)
# ---------------------------------------------------------------------------
FACE_TARGETS = [10000, 20000, 50000]   # decimation targets to characterise
PROD_FACES   = 20000                   # production decimation target
PRE_DECIMATE_THRESHOLD = 3_000_000     # SPHARM_main.py pre-decimates above this


def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "_targets.R").exists():
            return p
    return start.parents[2]


THIS_DIR  = Path(__file__).resolve().parent
PROJ_ROOT = find_project_root(THIS_DIR)
MESH_DIR  = PROJ_ROOT / "analysis" / "data" / "3D_models_cores"
MORPH_CSV = PROJ_ROOT / "analysis" / "data" / "derived_data" / "SPHARM_morphology.csv"
OUT_DIR   = THIS_DIR
FIG_DIR   = OUT_DIR / "figures"


def assemblage_of(idstr: str) -> str:
    s = str(idstr).strip()
    return "IM" if s.startswith("IM_") else ("SDG" if s.startswith("SDG")
            else ("EXP" if s.startswith("EXP") else "OTHER"))


def stl_header_info(path: Path):
    """Return (format, n_faces_header, file_bytes). ASCII face count is not read
    from the header (it would require a full parse); we rely on the committed
    n_faces_original for counts and use the header only for the format flag."""
    size = path.stat().st_size
    with open(path, "rb") as f:
        head = f.read(80)
        is_ascii = head.lstrip().startswith(b"solid")
        n_hdr = None
        if not is_ascii:
            raw = f.read(4)
            if len(raw) == 4:
                n_hdr = struct.unpack("<I", raw)[0]
    return ("ascii" if is_ascii else "binary"), n_hdr, size


def main() -> None:
    print("=" * 70)
    print("Mesh inventory (decimation context)")
    print(f"  project root : {PROJ_ROOT}")
    print(f"  meshes       : {MESH_DIR.relative_to(PROJ_ROOT)}")
    print(f"  face targets : {FACE_TARGETS}  (production = {PROD_FACES})")
    print("=" * 70)

    # Authoritative raw face counts come from the committed morphology cache
    # (n_faces_original, written by SPHARM_main.py as the pre-decimation count).
    morph = None
    if MORPH_CSV.exists():
        morph = pd.read_csv(MORPH_CSV, usecols=lambda c: c in ("ID", "n_faces_original"))
        morph["ID"] = morph["ID"].astype(str).str.strip()

    rows = []
    # Case-insensitive, de-duplicated listing (a plain glob("*.stl")+glob("*.STL")
    # double-counts on case-insensitive filesystems such as Windows/macOS).
    stl_files = sorted({p for p in MESH_DIR.iterdir()
                        if p.is_file() and p.suffix.lower() == ".stl"},
                       key=lambda p: p.name)
    if not stl_files:
        print(f"  [warn] no STL files found under {MESH_DIR}.")
    for p in stl_files:
        sid = p.stem.strip()
        fmt, n_hdr, size = stl_header_info(p)
        rows.append({"ID": sid, "group": assemblage_of(sid), "stl_format": fmt,
                     "n_faces_header": n_hdr, "file_mb": round(size / 1e6, 2)})
    inv = pd.DataFrame(rows)

    # Merge the authoritative raw face count.
    if morph is not None and not inv.empty:
        inv = inv.merge(morph, on="ID", how="left")
        inv["n_faces_raw"] = inv["n_faces_original"].fillna(inv["n_faces_header"])
    elif not inv.empty:
        inv["n_faces_raw"] = inv["n_faces_header"]

    # If meshes are absent (OSF-only), fall back to the cache alone so the report
    # still works from committed data.
    if inv.empty and morph is not None:
        inv = morph.assign(group=morph["ID"].map(assemblage_of),
                           stl_format="(mesh not present)",
                           n_faces_header=np.nan, file_mb=np.nan,
                           n_faces_raw=morph["n_faces_original"])

    if inv.empty:
        raise FileNotFoundError("No mesh files and no morphology cache to inventory.")

    inv = inv.sort_values(["group", "ID"]).reset_index(drop=True)
    inv.to_csv(OUT_DIR / "mesh_inventory.csv", index=False)

    # Per-group x target: how many specimens are actually decimated (raw > target)
    # vs returned unchanged (raw <= target).
    groups = ["IM", "SDG", "EXP"]
    summ_rows = []
    for grp in groups:
        g = inv[inv["group"] == grp]
        if g.empty:
            continue
        fr = g["n_faces_raw"].dropna()
        for T in FACE_TARGETS:
            dec = int((fr > T).sum())
            unch = int((fr <= T).sum())
            summ_rows.append({
                "group": grp, "face_target": T,
                "n_specimens": int(len(g)),
                "n_decimated": dec, "n_unchanged": unch,
                "pct_decimated": round(100.0 * dec / len(fr), 1) if len(fr) else np.nan,
                "raw_faces_min": int(fr.min()) if len(fr) else None,
                "raw_faces_median": int(fr.median()) if len(fr) else None,
                "raw_faces_max": int(fr.max()) if len(fr) else None,
                "n_pre_decimated_gt3M": int((fr > PRE_DECIMATE_THRESHOLD).sum()),
            })
    summ = pd.DataFrame(summ_rows)
    summ.to_csv(OUT_DIR / "mesh_inventory_summary.csv", index=False)

    print("\nRaw-mesh inventory by group:")
    for grp in groups:
        g = inv[inv["group"] == grp]
        if g.empty:
            continue
        fr = g["n_faces_raw"].dropna()
        fmts = g["stl_format"].value_counts().to_dict()
        print(f"  {grp}: n={len(g)}  raw faces min/median/max = "
              f"{int(fr.min()):,}/{int(fr.median()):,}/{int(fr.max()):,}  formats={fmts}")
    print("\nDecimation reach (specimens actually reduced at each target):")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.to_string(index=False))

    make_figure(inv, groups)

    print("\nWrote mesh_inventory.csv, mesh_inventory_summary.csv")
    print("Done.")


def make_figure(inv: pd.DataFrame, groups) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  [info] matplotlib unavailable ({e}); skipping inventory figure.")
        return False
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    palette = {"IM": "#5C7F71", "SDG": "#4A6E8A", "EXP": "#BA8530"}
    fig, ax = plt.subplots(figsize=(9, 4.6))
    for grp in groups:
        g = inv[inv["group"] == grp]
        fr = g["n_faces_raw"].dropna()
        if fr.empty:
            continue
        y = np.random.default_rng(0).normal(loc={"IM": 0, "SDG": 1, "EXP": 2}[grp],
                                             scale=0.06, size=len(fr))
        ax.scatter(fr, y, s=14, alpha=0.7, color=palette.get(grp, "grey"), label=grp)
    for T in FACE_TARGETS:
        ax.axvline(T, ls="--" if T != PROD_FACES else "-",
                   lw=1.0 if T == PROD_FACES else 0.7, color="grey")
        ax.text(T, 2.55, f"{T//1000}k", ha="center", va="bottom", fontsize=8, color="grey")
    ax.set_xscale("log")
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["IM", "SDG", "EXP"])
    ax.set_xlabel("Raw mesh face count (log scale)")
    ax.legend(title="group", frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIG_DIR / "fig_S_mesh_inventory.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(PROJ_ROOT)}")
    return True


if __name__ == "__main__":
    main()
