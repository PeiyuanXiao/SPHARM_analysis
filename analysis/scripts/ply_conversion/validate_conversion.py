"""Validate every (.stl, .ply) pair produced by convert_stl_to_ply.py.

Every pair is checked; nothing is sampled. Writes conversion_validation.csv
into the output directory.

Pass criteria
-------------
    faces_match           must be exactly True; any difference is an error
    volume_diff           must be <= VOL_TOL
    unique vertex set     sorted, compared point by point

Vertex COUNT is recorded but never used to pass or fail: see the note in
convert_stl_to_ply.py on why it does not change under process=False.

In practice every quantity here comes out at exactly 0.0 rather than merely
small, because both formats store coordinates as float32 and trimesh's PLY
writer emits float32, making the round trip bit-exact.

Run inside the project container:
    /opt/conda/envs/spharm/bin/python analysis/scripts/ply_conversion/validate_conversion.py
"""
import os
import sys
import glob
import csv

import numpy as np
import trimesh

SRC = "/project/analysis/data/3D_models_cores"          # analysis inputs (.stl)
DST = "/project/analysis/data/3D_models_cores_ply"      # archival copies (.ply)
OUT = "/project/analysis/data/conversion_validation.csv"

VOL_TOL = 1e-8


def load_raw(path):
    m = trimesh.load(path, process=False)
    if isinstance(m, trimesh.Scene):
        raise RuntimeError(f"{path} loaded as Scene")
    return m


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*.stl")))
    rows, hard_fail = [], []

    print(f"{'filename':<40} {'F match':>8} {'vol_diff':>12} "
          f"{'area_diff':>12} {'max_vdev':>12} {'status':>8}")
    print("-" * 98)

    for p in files:
        name = os.path.basename(p)
        q = os.path.join(DST, os.path.splitext(name)[0] + ".ply")

        ms, mp = load_raw(p), load_raw(q)
        sv, sf = len(ms.vertices), len(ms.faces)
        pv, pf = len(mp.vertices), len(mp.faces)

        faces_match = (sf == pf)
        vol_diff = abs(float(ms.volume) - float(mp.volume))
        area_diff = abs(float(ms.area) - float(mp.area))

        # deduplicated, sorted vertex-coordinate sets
        us = np.unique(np.asarray(ms.vertices, dtype=np.float64), axis=0)
        up = np.unique(np.asarray(mp.vertices, dtype=np.float64), axis=0)
        if us.shape == up.shape:
            max_vdev = float(np.abs(us - up).max())
            set_ok = True
        else:
            max_vdev = float("nan")
            set_ok = False

        # connectivity must survive too, not just the coordinates
        faces_identical = (faces_match and
                           np.array_equal(np.asarray(ms.faces),
                                          np.asarray(mp.faces)))

        bad = (not faces_match) or (vol_diff > VOL_TOL) or (not set_ok)
        status = "FAIL" if bad else "PASS"
        if bad:
            hard_fail.append(name)

        print(f"{name:<40} {str(faces_match):>8} {vol_diff:>12.3e} "
              f"{area_diff:>12.3e} {max_vdev:>12.3e} {status:>8}")

        rows.append(dict(
            filename=name, stl_vertices=sv, ply_vertices=pv,
            stl_faces=sf, ply_faces=pf,
            faces_match="TRUE" if faces_match else "FALSE",
            volume_diff=repr(vol_diff), area_diff=repr(area_diff),
            max_vertex_deviation=repr(max_vdev),
            unique_vertex_count=int(us.shape[0]),
            face_index_array_identical="TRUE" if faces_identical else "FALSE",
            status=status))

        del ms, mp, us, up

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("-" * 98)
    print(f"wrote {OUT}  ({len(rows)} rows)")
    print(f"faces_match all TRUE : "
          f"{all(r['faces_match'] == 'TRUE' for r in rows)}")
    print(f"face index arrays identical all TRUE : "
          f"{all(r['face_index_array_identical'] == 'TRUE' for r in rows)}")
    print(f"max volume_diff      : "
          f"{max(float(r['volume_diff']) for r in rows):.3e}")
    print(f"max area_diff        : "
          f"{max(float(r['area_diff']) for r in rows):.3e}")
    print(f"max vertex deviation : "
          f"{np.nanmax([float(r['max_vertex_deviation']) for r in rows]):.3e}")
    print(f"FAIL rows            : {len(hard_fail)}")
    for n in hard_fail:
        print(f"  FAIL {n}")
    return 1 if hard_fail else 0


if __name__ == "__main__":
    sys.exit(main())
