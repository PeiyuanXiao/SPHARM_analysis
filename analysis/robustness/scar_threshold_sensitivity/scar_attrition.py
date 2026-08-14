"""scar_attrition.py"""

from __future__ import annotations

import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# PARAMETERS  (keep THRESHOLDS in sync with 01_*.py and 02_*.R)
# ---------------------------------------------------------------------------
THRESHOLDS = [0.0, 5.0, 10.0]
MIN_VIABLE_SCARS = 3

APPLIED_GROUPS = ["EXP", "SDG"]


# ---------------------------------------------------------------------------
# Paths (discover project root via _targets.R, like sweep_spharm_bandwidth.py)
# ---------------------------------------------------------------------------
def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "_targets.R").exists():
            return p
    return start.parents[2]


THIS_DIR  = Path(__file__).resolve().parent
PROJ_ROOT = find_project_root(THIS_DIR)
RAW_XLSX  = PROJ_ROOT / "analysis" / "data" / "raw_data" / "Scar_orientation_data.xlsx"
OUT_DIR   = THIS_DIR
FIG_DIR   = OUT_DIR / "figures"

SHEETS = {0: "IM", 1: "SDG", 2: "EXP"}


def assemblage_of(idstr: str) -> str:
    s = str(idstr).strip()
    if s.startswith("IM_"):
        return "IM"
    if s.startswith("SDG"):
        return "SDG"
    if s.startswith("EXP"):
        return "EXP"
    return "OTHER"


# ---------------------------------------------------------------------------
# Load raw scars and compute lengths
# ---------------------------------------------------------------------------
def load_scars() -> pd.DataFrame:
    if not RAW_XLSX.exists():
        raise FileNotFoundError(f"Raw scar data not found: {RAW_XLSX}")
    xl = pd.ExcelFile(RAW_XLSX)
    frames = []
    for idx in SHEETS:
        df = xl.parse(xl.sheet_names[idx])
        if df.empty:
            continue
        frames.append(df)
    raw = pd.concat(frames, ignore_index=True)

    raw["ID"] = raw["ID"].astype(str).str.strip()

    need = {"Start_X", "Start_Y", "Start_Z", "End_X", "End_Y", "End_Z"}
    missing = need - set(raw.columns)
    if missing:
        raise ValueError(f"Workbook missing endpoint columns: {missing}")

    raw["length_mm"] = np.sqrt(
        (raw["End_X"] - raw["Start_X"]) ** 2
        + (raw["End_Y"] - raw["Start_Y"]) ** 2
        + (raw["End_Z"] - raw["Start_Z"]) ** 2
    )
    raw["group"] = raw["ID"].map(assemblage_of)
    if "Typology" not in raw.columns:
        raw["Typology"] = ""
    raw["Typology"] = raw["Typology"].fillna("").astype(str).str.strip()
    return raw


# ---------------------------------------------------------------------------
# Per-specimen and per-assemblage attrition tables
# ---------------------------------------------------------------------------
def per_specimen_table(raw: pd.DataFrame) -> pd.DataFrame:
    rows = []
    typ = (raw.groupby("ID")["Typology"]
              .agg(lambda s: next((x for x in s if x), ""))
              .to_dict())
    for (gid, idv), g in raw.groupby(["group", "ID"]):
        n_total = len(g)
        for T in THRESHOLDS:
            kept = g.loc[g["length_mm"] > T]
            rows.append({
                "ID": idv, "group": gid, "Typology": typ.get(idv, ""),
                "threshold_mm": T,
                "n_total": n_total,
                "n_retained": int(len(kept)),
                "pct_retained": round(100.0 * len(kept) / n_total, 2) if n_total else np.nan,
                "L_min_retained": round(float(kept["length_mm"].min()), 3) if len(kept) else np.nan,
                "L_max_retained": round(float(kept["length_mm"].max()), 3) if len(kept) else np.nan,
            })
    out = pd.DataFrame(rows).sort_values(["group", "ID", "threshold_mm"])
    return out.reset_index(drop=True)


def summary_table(per_spec: pd.DataFrame, groups) -> pd.DataFrame:
    rows = []
    for grp in groups:
        gs = per_spec[per_spec["group"] == grp]
        n_spec = gs["ID"].nunique()
        for T in THRESHOLDS:
            sub = gs[gs["threshold_mm"] == T]
            tot = int(sub["n_total"].sum())
            ret = int(sub["n_retained"].sum())
            rc = sub["n_retained"]
            rows.append({
                "group": grp,
                "threshold_mm": T,
                "n_specimens": int(n_spec),
                "scars_total": tot,
                "scars_retained": ret,
                "scars_dropped": tot - ret,
                "pct_retained": round(100.0 * ret / tot, 2) if tot else np.nan,
                "per_spec_min": int(rc.min()),
                "per_spec_median": float(np.median(rc)),
                "per_spec_max": int(rc.max()),
                f"n_spec_below_{MIN_VIABLE_SCARS}": int((rc < MIN_VIABLE_SCARS).sum()),
                "n_spec_zero": int((rc == 0).sum()),
            })
    return pd.DataFrame(rows)


def im_diagnostic(per_spec: pd.DataFrame) -> pd.DataFrame:
    """Document WHY the threshold is not applied to the synthetic ideal cores."""
    im = per_spec[per_spec["group"] == "IM"].copy()
    wide = im.pivot_table(index=["ID", "n_total"], columns="threshold_mm",
                          values="n_retained", aggfunc="first").reset_index()
    wide.columns = ["ID", "n_total"] + [f"n_retained_gt{int(c)}mm" for c in THRESHOLDS]
    lmin = im.groupby("ID")["L_min_retained"].min()
    raw_len = per_spec
    wide["note"] = np.where(
        wide.filter(like="n_retained_gt5mm").iloc[:, 0] == 0,
        "ERASED at >5mm (synthetic uniform length)", "")
    return wide


# ---------------------------------------------------------------------------
# Figure (best-effort: needs matplotlib; the canonical SI figure is built by 02_*.R)
# ---------------------------------------------------------------------------
def make_figure(per_spec: pd.DataFrame) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  [info] matplotlib unavailable ({e}); skipping local preview figure."
              "\n         The canonical scar-count figure is generated by 02_*.R (ggplot).")
        return False

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    palette = {0.0: "#4A6E8A", 5.0: "#BA8530", 10.0: "#802520"}
    fig, axes = plt.subplots(1, len(APPLIED_GROUPS), figsize=(11, 4.6), squeeze=False)
    for j, grp in enumerate(APPLIED_GROUPS):
        ax = axes[0][j]
        gs = per_spec[per_spec["group"] == grp]
        base = (gs[gs["threshold_mm"] == 0.0]
                .sort_values("n_retained", ascending=False)["ID"].tolist())
        rank = {idv: r for r, idv in enumerate(base)}
        for T in THRESHOLDS:
            sub = gs[gs["threshold_mm"] == T].copy()
            sub["rank"] = sub["ID"].map(rank)
            sub = sub.sort_values("rank")
            ax.plot(sub["rank"], sub["n_retained"], marker="o", ms=3, lw=1.2,
                    color=palette.get(T, "grey"),
                    label=("all scars (~2 mm)" if T == 0 else f"> {int(T)} mm"))
        ax.axhline(MIN_VIABLE_SCARS, ls="--", lw=0.6, color="grey")
        ax.set_title(f"{grp} cores  (n = {gs['ID'].nunique()})")
        ax.set_xlabel("specimen (ranked by baseline scar count)")
        if j == 0:
            ax.set_ylabel("scars retained")
        ax.legend(title="minimum size", fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIG_DIR / "fig_S_threshold_attrition.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(PROJ_ROOT)}")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 70)
    print("Scar attrition report (minimum-size thresholds)")
    print(f"  project root : {PROJ_ROOT}")
    print(f"  raw data     : {RAW_XLSX.relative_to(PROJ_ROOT)}")
    print(f"  thresholds   : {THRESHOLDS} mm  (keep length > T; T=0 = all recorded)")
    print("=" * 70)

    raw = load_scars()
    n_by_grp = raw.groupby("group")["ID"].nunique().to_dict()
    print(f"Loaded {len(raw)} scars across {raw['ID'].nunique()} specimens "
          f"(EXP={n_by_grp.get('EXP',0)}, SDG={n_by_grp.get('SDG',0)}, "
          f"IM={n_by_grp.get('IM',0)}).\n")

    per_spec = per_specimen_table(raw)
    summ     = summary_table(per_spec, APPLIED_GROUPS)
    im_diag  = im_diagnostic(per_spec)

    per_spec.to_csv(OUT_DIR / "scar_attrition_by_specimen.csv", index=False)
    summ.to_csv(OUT_DIR / "scar_attrition_summary.csv", index=False)
    im_diag.to_csv(OUT_DIR / "im_threshold_diagnostic.csv", index=False)
    print("Wrote scar_attrition_by_specimen.csv, scar_attrition_summary.csv, "
          "im_threshold_diagnostic.csv")

    print("\nRetention summary (EXP, SDG):")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ.to_string(index=False))
    print("\nIdeal-core diagnostic (threshold NOT applied; held at production):")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(im_diag.to_string(index=False))

    make_figure(per_spec)

    print("\nDone.")


if __name__ == "__main__":
    main()
