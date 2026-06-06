"""
01_sweep_mesh_spharm.py
=======================
Mesh-PREPROCESSING sensitivity sweep for the M-SPHARM (morphology) pipeline.
NEW, self-contained add-on for the paper's Supplementary Information. It does NOT
modify the main pipeline, the cached results, or the manuscript, and it RE-USES the
existing mesh / spherical-harmonic functions unchanged (open3d quadric decimation,
trimesh Laplacian smoothing, and the project's `mesh_processing`, `pca_align`,
`spherical_harmonics`, `power_spectrum` modules). The orchestration below mirrors
`SPHARM_main.process_single_mesh` verbatim; the ONLY things that vary across the
sweep are the two preprocessing parameters under test:

    (1) the decimation TARGET FACE COUNT  (production = 20000)
    (2) the Laplacian SMOOTHING iterations (production = 3)

The decimation ALGORITHM is held fixed at quadric edge collapse
(`open3d.simplify_quadric_decimation`); every other step — pre-decimation policy,
trimesh cleaning, `mesh_processing.normalize_mesh`, `pca_align.robust_pca_alignment`
(deterministic, area-weighted PCA + sign convention), the GRID_SIZE = 256 spherical
interpolation, the l_max = 20 expansion, 4pi/csphase normalisation, and all seeds —
is identical to the main pipeline. The scar (SP-SPHARM) side does not derive from
these meshes (scar vectors come from digitised endpoint coordinates), so it is
unaffected and is reused unchanged downstream.

WHAT IT DOES, for each (face_target, smooth_iters) in SETTINGS:
    each core mesh .stl  ->  process_mesh(target_faces, smooth_iters)
        -> spectra/SPHARM_morphology_f{face}_s{smooth}.csv
           (ID, n_faces_original, n_faces_decimated, power_l0 .. power_l20)

TWO BUILT-IN CHECKS at the production setting (20000, 3):
  * ANCHOR — the recomputed M-SPHARM power spectra must reproduce the cached
    production file derived_data/SPHARM_morphology.csv (max|diff| recorded in
    sweep_manifest.csv; flagged if > 1e-6, e.g. a BLAS / open3d / library mismatch,
    pinned in analysis/scripts/environment.yml).
  * DETERMINISM — re-running the same input through the fixed decimation twice must
    return an identical mesh AND identical spectra. Historically, manual/inconsistent
    decimation was the source of non-reproducible morphology; this check confirms the
    fixed-algorithm pipeline is now deterministic. Result in determinism_check.csv.
    A mismatch here is itself an important finding (decimation still not deterministic).

HOW TO RUN (canonical environment — conda `spharm`, same as the main pipeline):
    python analysis/mesh_decimation_sensitivity/01_sweep_mesh_spharm.py
    # or, mirroring _targets.R: PYTHONPATH=analysis/scripts python .../01_sweep_mesh_spharm.py

Then run 02_mesh_sensitivity_stats.R to evaluate stability of the downstream
conclusions across settings.
"""

import gc
import os
import struct
import sys
import tempfile
import time
import platform
from pathlib import Path

# ---------------------------------------------------------------------------
# PARAMETERS  (keep in sync with 00_mesh_inventory.py and 02_*.R)
# ---------------------------------------------------------------------------
FACE_TARGETS = [10000, 20000, 50000]   # decimation target face counts (production = 20000)
SMOOTH_ITERS = [0, 3, 6]               # Laplacian smoothing iterations (production = 3)
PROD_FACES   = 20000
PROD_SMOOTH  = 3

# Two 1-D sweeps through the production setting (vary faces at production smoothing;
# vary smoothing at production faces). Set FULL_CROSS=True to also run the full grid.
FULL_CROSS = False

# Fixed pipeline settings — IDENTICAL to SPHARM_main.py.
GRID_SIZE = 256                        # spherical-interpolation grid (MUST equal 256 to match production)
LMAX      = 20
PRE_DECIMATE_THRESHOLD = 3_000_000     # pre-decimate above this face count
PRE_DECIMATE_TARGET    = 500_000

# Cost controls. SPECIMEN_LIMIT=None processes all meshes; set an int for a smoke test.
# N_DETERMINISM=None checks determinism on every specimen at the production setting.
SPECIMEN_LIMIT = None
N_DETERMINISM  = None


def build_settings():
    s = [(f, PROD_SMOOTH) for f in FACE_TARGETS] + [(PROD_FACES, it) for it in SMOOTH_ITERS]
    if FULL_CROSS:
        s += [(f, it) for f in FACE_TARGETS for it in SMOOTH_ITERS]
    seen, out = set(), []
    for st in s:
        if st not in seen:
            seen.add(st); out.append(st)
    return out


# ---------------------------------------------------------------------------
# Paths + module imports (mirror _targets.R / SPHARM_main.py)
# ---------------------------------------------------------------------------
def find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "_targets.R").exists():
            return p
    return start.parents[2]


THIS_DIR    = Path(__file__).resolve().parent
PROJ_ROOT   = find_project_root(THIS_DIR)
SCRIPTS_DIR = PROJ_ROOT / "analysis" / "scripts"
MESH_DIR    = PROJ_ROOT / "analysis" / "data" / "3D_models_cores"
DERIVED_DIR = PROJ_ROOT / "analysis" / "data" / "derived_data"
OUT_DIR     = THIS_DIR
SPECTRA_DIR = OUT_DIR / "spectra"
CACHE_CSV   = DERIVED_DIR / "SPHARM_morphology.csv"   # cached production (20000,3) file

sys.path.insert(0, str(SCRIPTS_DIR))

import numpy as np
import pandas as pd
import pyshtools as pysh
import trimesh
from trimesh.smoothing import filter_laplacian

from SPHARM_modules import mesh_processing, pca_align, spherical_harmonics
from SPHARM_modules.power_spectrum import compute_power_spectrum

np.random.seed(42)
POWER_COLS = [f"power_l{l}" for l in range(LMAX + 1)]


# ---------------------------------------------------------------------------
# Core: one mesh -> M-SPHARM, parameterised by (target_faces, smooth_iters).
# This mirrors SPHARM_main.process_single_mesh exactly; only TARGET_FACES and the
# Laplacian iteration count are turned into parameters, and the decimated mesh is
# returned alongside the feature row for the determinism check.
# ---------------------------------------------------------------------------
def process_mesh(stl_path, target_faces: int, smooth_iters: int):
    import open3d as o3d

    specimen_id = os.path.splitext(os.path.basename(stl_path))[0]
    tmp_path = None

    with open(stl_path, "rb") as f:
        header_bytes = f.read(80)
        is_ascii = header_bytes.lstrip().startswith(b"solid")
        if not is_ascii:
            n_raw = struct.unpack("<I", f.read(4))[0]
        else:
            n_raw = None

    if is_ascii:
        o3d_mesh = o3d.io.read_triangle_mesh(stl_path)
        n_raw = len(o3d_mesh.triangles)
        if n_raw > PRE_DECIMATE_THRESHOLD:
            o3d_mesh = o3d_mesh.simplify_quadric_decimation(PRE_DECIMATE_TARGET)
            o3d_mesh.compute_vertex_normals()
            tmp_path = tempfile.mktemp(suffix=".stl")
            o3d.io.write_triangle_mesh(tmp_path, o3d_mesh)
            del o3d_mesh; gc.collect()
            load_path = tmp_path
        else:
            load_path = None
            o3d_mesh_hold = o3d_mesh
    else:
        if n_raw > PRE_DECIMATE_THRESHOLD:
            tmp_path = tempfile.mktemp(suffix=".stl")
            step     = max(1, n_raw // PRE_DECIMATE_TARGET)
            keep_set = set(range(0, n_raw, step))
            n_keep   = len(keep_set)
            with open(stl_path, "rb") as fin, open(tmp_path, "wb") as fout:
                header = fin.read(80); fin.read(4)
                fout.write(header); fout.write(struct.pack("<I", n_keep))
                for i in range(n_raw):
                    face_data = fin.read(50)
                    if i in keep_set:
                        fout.write(face_data)
            gc.collect()
            load_path = tmp_path
        else:
            load_path = stl_path

    try:
        if is_ascii and n_raw <= PRE_DECIMATE_THRESHOLD:
            o3d_mesh_final = o3d_mesh_hold
        else:
            o3d_mesh_final = o3d.io.read_triangle_mesh(load_path)

        # (1) Decimation — quadric edge collapse to the swept target face count.
        o3d_mesh_final = o3d_mesh_final.simplify_quadric_decimation(target_faces)
        decimated_vertices = np.asarray(o3d_mesh_final.vertices)
        decimated_faces    = np.asarray(o3d_mesh_final.triangles)
        if len(decimated_faces) == 0:
            raise ValueError("open3d.simplify_quadric_decimation returned empty mesh")

        # (2) Clean + Laplacian smoothing (iterations swept; 0 = no smoothing).
        valid_mask = np.all(decimated_faces < len(decimated_vertices), axis=1)
        decimated_faces = decimated_faces[valid_mask]
        mesh = trimesh.Trimesh(vertices=decimated_vertices,
                               faces=decimated_faces, process=True)
        mesh.remove_unreferenced_vertices()
        if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
            raise ValueError("Mesh empty after cleaning unreferenced vertices")
        if smooth_iters > 0:
            filter_laplacian(mesh, iterations=smooth_iters, volume_constraint=False)
        decimated_vertices = mesh.vertices
        decimated_faces    = mesh.faces

        # (3) Normalise + PCA alignment (deterministic).
        normalized_vertices = mesh_processing.normalize_mesh(decimated_vertices, decimated_faces)
        aligned_vertices, _ = pca_align.robust_pca_alignment(
            normalized_vertices, faces=decimated_faces, enforce_direction=True)

        # (4) Spherical interpolation + SH expansion + power spectrum.
        spherical_coords = spherical_harmonics.cartesian_to_spherical(aligned_vertices)
        R, theta, phi    = spherical_coords.T
        grid_r           = spherical_harmonics.spherical_interpolate(R, theta, phi, GRID_SIZE)
        clm    = spherical_harmonics.compute_spherical_harmonics(
                     grid_r, normalization_method="zero-component")
        clm_sh = pysh.SHCoeffs.from_array(clm, normalization="4pi", csphase=1, lmax=LMAX).pad(lmax=LMAX)
        norm_power = compute_power_spectrum(clm_sh, lmax=LMAX)["norm_power"]

        row = {"ID": specimen_id,
               "n_faces_original":  int(n_raw),
               "n_faces_decimated": int(len(decimated_faces))}
        for l, p in enumerate(norm_power):
            row[f"power_l{l}"] = float(p)
        return row, np.asarray(decimated_vertices), np.asarray(decimated_faces)

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


# ---------------------------------------------------------------------------
def list_meshes():
    files = sorted({p for p in MESH_DIR.iterdir()
                    if p.is_file() and p.suffix.lower() == ".stl"}, key=lambda p: p.name)
    if SPECIMEN_LIMIT:
        files = files[:SPECIMEN_LIMIT]
    return files


def spectra_for_setting(stl_files, target_faces, smooth_iters):
    rows = []
    for i, stl in enumerate(stl_files):
        sid = stl.stem
        try:
            row, _, _ = process_mesh(str(stl), target_faces, smooth_iters)
            rows.append(row)
        except Exception as e:                       # match batch_process robustness
            print(f"    [warn] {sid}: {e}")
            rows.append({"ID": sid})
        if (i + 1) % 20 == 0:
            print(f"      {i + 1}/{len(stl_files)} meshes")
    return pd.DataFrame(rows)


def max_abs_diff_vs_cache(df_s: pd.DataFrame) -> float:
    if not CACHE_CSV.exists():
        print(f"    [info] cache not found ({CACHE_CSV.name}); skipping anchor check")
        return float("nan")
    cache = pd.read_csv(CACHE_CSV)
    cache["ID"] = cache["ID"].astype(str).str.strip()
    cols = [c for c in POWER_COLS if c in cache.columns and c in df_s.columns]
    a = df_s.assign(ID=df_s["ID"].astype(str).str.strip()).set_index("ID")[cols].sort_index()
    b = cache.set_index("ID")[cols].reindex(a.index)
    return float(np.nanmax(np.abs(a.values - b.values)))


def determinism_check(stl_files):
    """Run the production setting twice per specimen; compare meshes and spectra."""
    sub = stl_files if N_DETERMINISM is None else stl_files[:N_DETERMINISM]
    print(f"\n[determinism] re-running ({PROD_FACES}, {PROD_SMOOTH}) twice on "
          f"{len(sub)} specimens...")
    out = []
    for stl in sub:
        sid = stl.stem
        try:
            r1, v1, f1 = process_mesh(str(stl), PROD_FACES, PROD_SMOOTH)
            r2, v2, f2 = process_mesh(str(stl), PROD_FACES, PROD_SMOOTH)
            mesh_identical = (v1.shape == v2.shape and f1.shape == f2.shape
                              and np.array_equal(v1, v2) and np.array_equal(f1, f2))
            p1 = np.array([r1.get(c, np.nan) for c in POWER_COLS])
            p2 = np.array([r2.get(c, np.nan) for c in POWER_COLS])
            pdiff = float(np.nanmax(np.abs(p1 - p2)))
            out.append({"ID": sid, "mesh_identical": bool(mesh_identical),
                        "n_faces_decimated": r1.get("n_faces_decimated"),
                        "max_power_diff": pdiff})
        except Exception as e:
            out.append({"ID": sid, "mesh_identical": None,
                        "n_faces_decimated": None, "max_power_diff": float("nan"),
                        "error": str(e)})
    df = pd.DataFrame(out)
    df.to_csv(OUT_DIR / "determinism_check.csv", index=False)
    n_ok = int(df["mesh_identical"].fillna(False).sum())
    max_pdiff = float(np.nanmax(df["max_power_diff"])) if len(df) else float("nan")
    deterministic = (n_ok == len(df)) and (max_pdiff <= 1e-9)
    print(f"  meshes identical on re-run: {n_ok}/{len(df)};  max spectra diff = {max_pdiff:.2e}")
    print(f"  => {'DETERMINISTIC (fixed-algorithm decimation reproduces exactly).' if deterministic else '*** NON-DETERMINISTIC — investigate (this is a reportable finding).'}")
    return {"n_identical": n_ok, "n_total": int(len(df)),
            "max_power_diff": max_pdiff, "deterministic": bool(deterministic)}


# ---------------------------------------------------------------------------
def main() -> None:
    settings = build_settings()
    print("=" * 70)
    print("M-SPHARM mesh-preprocessing sensitivity sweep")
    print(f"  project root : {PROJ_ROOT}")
    print(f"  meshes       : {MESH_DIR.relative_to(PROJ_ROOT)}")
    print(f"  face targets : {FACE_TARGETS} (prod {PROD_FACES})")
    print(f"  smooth iters : {SMOOTH_ITERS} (prod {PROD_SMOOTH})")
    print(f"  settings     : {settings}")
    print(f"  grid/lmax    : GRID_SIZE={GRID_SIZE} / lmax={LMAX}")
    print("=" * 70)

    if not MESH_DIR.exists():
        raise FileNotFoundError(
            f"Mesh directory not found: {MESH_DIR}\n"
            "Download 3D_models_cores from OSF (see README / Data section) first.")
    SPECTRA_DIR.mkdir(parents=True, exist_ok=True)
    stl_files = list_meshes()
    print(f"Found {len(stl_files)} meshes.\n")

    manifest = []
    for (face, smooth) in settings:
        is_prod = (face == PROD_FACES and smooth == PROD_SMOOTH)
        print(f"[faces={face}, smooth={smooth}]{'  (PRODUCTION)' if is_prod else ''}")
        t0 = time.time()
        df_s = spectra_for_setting(stl_files, face, smooth)
        elapsed = time.time() - t0

        out_csv = SPECTRA_DIR / f"SPHARM_morphology_f{face}_s{smooth}.csv"
        df_s.to_csv(out_csv, index=False)

        diff, flag = float("nan"), ""
        if is_prod:
            diff = max_abs_diff_vs_cache(df_s)
            if not np.isnan(diff):
                ok = diff <= 1e-6
                flag = "OK" if ok else "*** MISMATCH ***"
                print(f"    anchor vs cache: max|diff| = {diff:.3e}  [{flag}]")
                if not ok:
                    print("    >>> Recomputed production spectra DIFFER from the cached "
                          "SPHARM_morphology.csv. Check open3d/trimesh/pyshtools/BLAS "
                          "versions against analysis/scripts/environment.yml.")

        manifest.append({
            "face_target": face, "smooth_iters": smooth, "is_production": is_prod,
            "n_specimens": int(df_s["ID"].nunique()),
            "elapsed_s": round(elapsed, 1),
            "sec_per_specimen": round(elapsed / max(len(stl_files), 1), 2),
            "output_csv": str(out_csv.relative_to(PROJ_ROOT)),
            "max_abs_diff_cache": diff, "anchor_flag": flag,
        })
        print(f"    wrote {out_csv.relative_to(PROJ_ROOT)}  ({elapsed:.1f}s)\n")

    # Determinism check at the production setting.
    det = determinism_check(stl_files)

    man_df = pd.DataFrame(manifest)
    man_df.to_csv(OUT_DIR / "sweep_manifest.csv", index=False)
    pd.DataFrame([det]).to_csv(OUT_DIR / "determinism_summary.csv", index=False)
    print(f"\nWrote sweep_manifest.csv and determinism_summary.csv")

    # Versions.
    def ver(mod):
        try:
            return __import__(mod).__version__
        except Exception:
            return "not importable"
    print("\nDone. Next: Rscript analysis/mesh_decimation_sensitivity/02_mesh_sensitivity_stats.R")


if __name__ == "__main__":
    main()
