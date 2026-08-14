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
