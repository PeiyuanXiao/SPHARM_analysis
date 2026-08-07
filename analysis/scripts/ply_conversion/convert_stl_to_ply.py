"""Convert the archived core models from .stl to .ply.

Produced in response to reviewer 1, who asked for .ply rather than .stl as the
archival format.

WHY process=False
-----------------
trimesh's default load path merges duplicate vertices and drops degenerate
faces. That would make the exported .ply a *cleaned* mesh rather than a
faithful re-encoding of the file the analysis actually consumed. Archival
copies must not be silently modified, so geometry processing is disabled.

Note that with process=False the vertex count does NOT fall. An STL has no
shared vertex table (V = 3F, each triangle storing its three corners
independently); disabling the merge preserves that layout, so the .ply keeps
all 3F vertices. A faithful re-encoding and a smaller vertex count are
mutually exclusive here.

⚠ The .ply are an archival copy, NOT a drop-in pipeline input. Open3D's STL
reader silently merges a small number of vertices while its PLY reader does
not, and simplify_quadric_decimation is greedy and topology-sensitive, so
running the pipeline off the .ply shifts normalised power by ~1e-4 absolute
relative to the published values. SPHARM_main.py also filters on '.stl' and
parses binary STL headers directly. Archive the .ply alongside the .stl; do
not substitute it. See verify_equivalence.py for the measurements behind this.

Run inside the project container:
    /opt/conda/envs/spharm/bin/python analysis/scripts/ply_conversion/convert_stl_to_ply.py
"""
import os
import sys
import glob
import time

import trimesh

SRC = "/project/analysis/data/3D_models_cores"          # analysis inputs (.stl)
DST = "/project/analysis/data/3D_models_cores_ply"      # archival copies (.ply)


def load_raw(path):
    """Load a mesh with no geometry processing whatsoever."""
    m = trimesh.load(path, process=False)
    if isinstance(m, trimesh.Scene):
        raise RuntimeError(f"{path} loaded as Scene, not a single mesh")
    return m


def main():
    os.makedirs(DST, exist_ok=True)

    files = sorted(glob.glob(os.path.join(SRC, "*.stl")))
    print(f"source dir : {SRC}")
    print(f"output dir : {DST}")
    print(f"files found: {len(files)}\n")
    print(f"{'filename':<40} {'stl_V':>10} {'stl_F':>10} "
          f"{'ply_V':>10} {'ply_F':>10}  {'sec':>6}")
    print("-" * 96)

    ok, failed = 0, []
    for p in files:
        name = os.path.basename(p)
        out = os.path.join(DST, os.path.splitext(name)[0] + ".ply")
        t0 = time.time()
        try:
            m = load_raw(p)
            sv, sf = len(m.vertices), len(m.faces)
            m.export(out)
            del m

            r = load_raw(out)
            pv, pf = len(r.vertices), len(r.faces)
            del r

            print(f"{name:<40} {sv:>10,} {sf:>10,} {pv:>10,} {pf:>10,}  "
                  f"{time.time() - t0:>6.1f}")
            ok += 1
        except Exception as e:
            print(f"{name:<40} FAILED: {e}")
            failed.append((name, str(e)))

    print("-" * 96)
    print(f"converted: {ok} / {len(files)}   failed: {len(failed)}")
    for n, e in failed:
        print(f"  FAIL {n}: {e}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
