"""Two independent checks that the .stl -> .ply conversion is lossless.

A. Geometric equivalence. For every pair, realise the actual triangle corner
   coordinates (vertices[faces]) and compare bitwise. This tests coordinates
   and connectivity together, so it is a complete geometric identity proof.
   Expected result: bit-exact for every pair.

B. Full M-SPHARM re-run on a few specimens, from the .stl and from the .ply,
   each compared against the published SPHARM_morphology.csv. This mirrors
   SPHARM_main.process_single_mesh for the sub-threshold branch.

   ⚠ Check B is EXPECTED TO DIFFER, and the difference is not a conversion
   defect. Open3D's STL reader silently merges a handful of vertices while its
   PLY reader merges none (e.g. SDG_L2_489: 153,276 vs 153,282), and
   simplify_quadric_decimation is greedy and topology-sensitive, so the tiny
   input difference cascades into a slightly different 20,000-face mesh.
   Measured shift: ~1e-4 absolute / ~1e-3 relative in normalised power.

   The two control columns are what make this readable: the chain is exactly
   deterministic (0.0 across repeat runs) and reproduces the published CSV
   from the .stl to ~1e-16, so the .ply-vs-.stl gap is attributable to the
   reader difference alone and to nothing in this code.

   Conclusion: the .ply is a verified lossless re-encoding of the analysed
   .stl (check A), but must not be substituted as pipeline input (check B).

Run inside the project container:
    /opt/conda/envs/spharm/bin/python analysis/scripts/ply_conversion/verify_equivalence.py
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import trimesh
import open3d as o3d
import pyshtools as pysh
from trimesh.smoothing import filter_laplacian

sys.path.insert(0, "/project/analysis/scripts")
from SPHARM_modules import mesh_processing, pca_align, spherical_harmonics
from SPHARM_modules.power_spectrum import compute_power_spectrum

SRC = "/project/analysis/data/3D_models_cores"          # analysis inputs (.stl)
DST = "/project/analysis/data/3D_models_cores_ply"      # archival copies (.ply)
CSV = "/project/analysis/data/derived_data/SPHARM_morphology.csv"

# must match SPHARM_main.py
TARGET_FACES, GRID_SIZE, LMAX = 20000, 256, 20

SAMPLES = ["SDG_L2_489", "SDG_L3_1066", "SDG_L4_3283"]


def check_geometric_equivalence():
    print("=" * 78)
    print("A. GEOMETRIC EQUIVALENCE  (per-triangle corner coordinates, all pairs)")
    print("=" * 78)

    worst, n_bitexact = 0.0, 0
    files = sorted(glob.glob(os.path.join(SRC, "*.stl")))
    for p in files:
        q = os.path.join(DST, os.path.splitext(os.path.basename(p))[0] + ".ply")
        a = trimesh.load(p, process=False)
        b = trimesh.load(q, process=False)
        ta = np.asarray(a.vertices)[np.asarray(a.faces)]   # (F, 3, 3)
        tb = np.asarray(b.vertices)[np.asarray(b.faces)]
        if ta.shape != tb.shape:
            print(f"  SHAPE MISMATCH {os.path.basename(p)}: {ta.shape} vs {tb.shape}")
            worst = float("inf")
        else:
            worst = max(worst, float(np.abs(ta - tb).max()))
            if np.array_equal(ta, tb):
                n_bitexact += 1
        del a, b, ta, tb

    print(f"pairs compared              : {len(files)}")
    print(f"bit-exact triangle arrays   : {n_bitexact} / {len(files)}")
    print(f"max |coordinate difference| : {worst:.3e}")
    return worst


def spharm_chain(path):
    """Mirrors SPHARM_main.process_single_mesh, sub-threshold branch."""
    m = o3d.io.read_triangle_mesh(path)
    m = m.simplify_quadric_decimation(TARGET_FACES)
    v = np.asarray(m.vertices)
    f = np.asarray(m.triangles)
    f = f[np.all(f < len(v), axis=1)]

    mesh = trimesh.Trimesh(vertices=v, faces=f, process=True)
    mesh.remove_unreferenced_vertices()
    filter_laplacian(mesh, iterations=3, volume_constraint=False)

    nv = mesh_processing.normalize_mesh(mesh.vertices, mesh.faces)
    av, _ = pca_align.robust_pca_alignment(nv, faces=mesh.faces,
                                           enforce_direction=True)
    R, theta, phi = spherical_harmonics.cartesian_to_spherical(av).T
    grid_r = spherical_harmonics.spherical_interpolate(R, theta, phi, GRID_SIZE)
    clm = spherical_harmonics.compute_spherical_harmonics(
        grid_r, normalization_method="zero-component")
    clm_sh = pysh.SHCoeffs.from_array(clm, normalization="4pi",
                                      csphase=1, lmax=LMAX).pad(lmax=LMAX)
    return np.asarray(compute_power_spectrum(clm_sh, lmax=LMAX)["norm_power"],
                      dtype=np.float64)


def check_pipeline_reruns():
    df = pd.read_csv(CSV).set_index("ID")
    pcols = [f"power_l{l}" for l in range(LMAX + 1)]

    print()
    print("=" * 78)
    print(f"B. FULL M-SPHARM RE-RUN  ({len(SAMPLES)} specimens)")
    print("=" * 78)

    def report(tag, a, b):
        ad = np.abs(a - b)
        with np.errstate(divide="ignore", invalid="ignore"):
            rd = np.where(np.abs(b) > 0, ad / np.abs(b), 0.0)
        print(f"    {tag:<34} max|abs|={ad.max():.6e}   "
              f"max|rel|={np.nanmax(rd):.6e}")
        return ad.max()

    summary = []
    for sid in SAMPLES:
        print(f"\n--- {sid} ---")
        pub = df.loc[sid, pcols].to_numpy(dtype=np.float64)
        stl = spharm_chain(os.path.join(SRC, sid + ".stl"))
        stl2 = spharm_chain(os.path.join(SRC, sid + ".stl"))
        ply = spharm_chain(os.path.join(DST, sid + ".ply"))

        summary.append((
            sid,
            report("stl run #1 vs stl run #2 (determinism)", stl, stl2),
            report("stl re-run vs published CSV", stl, pub),
            report("ply re-run vs published CSV", ply, pub),
            report("ply re-run vs stl re-run", ply, stl),
        ))

    print("\n" + "=" * 78)
    print(f"{'specimen':<16} {'determinism':>13} {'stl-vs-pub':>13} "
          f"{'ply-vs-pub':>13} {'ply-vs-stl':>13}")
    for sid, a, b, c, d in summary:
        print(f"{sid:<16} {a:>13.3e} {b:>13.3e} {c:>13.3e} {d:>13.3e}")


def main():
    worst = check_geometric_equivalence()
    check_pipeline_reruns()

    print()
    print("environment")
    print(f"  python    {sys.version.split()[0]}")
    print(f"  numpy     {np.__version__}")
    print(f"  trimesh   {trimesh.__version__}")
    print(f"  open3d    {o3d.__version__}")
    print(f"  pyshtools {pysh.__version__}")
    print(f"  pandas    {pd.__version__}")

    # Only check A is a pass/fail gate; see the module docstring for why B is
    # expected to differ.
    return 0 if worst == 0.0 else 1


if __name__ == "__main__":
    sys.exit(main())
