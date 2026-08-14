"""mesh_resolution_check.py"""
import os
import struct
import sys

import numpy as np
import pandas as pd
import open3d as o3d
import trimesh
from trimesh.smoothing import filter_laplacian

MESH_DIR = "/project/analysis/data/3D_models_cores"
OUT_CSV = "/project/analysis/robustness/sdg_projection_explore/mesh_resolution.csv"

TARGET_FACES = 20000
SMOOTH_ITERS = 3


def raw_face_count(path):
    """Face count from the STL itself: binary header, or a full parse if ASCII."""
    with open(path, "rb") as f:
        header = f.read(80)
        if header.lstrip().startswith(b"solid"):
            return int(len(o3d.io.read_triangle_mesh(path).triangles)), "ascii"
        return int(struct.unpack("<I", f.read(4))[0]), "binary"


def main():
    files = sorted(
        f for f in os.listdir(MESH_DIR)
        if f.lower().endswith(".stl") and not f.startswith("IM_")
    )
    print(f"measuring {len(files)} meshes (IM_ excluded)", flush=True)

    rows = []
    for i, fn in enumerate(files, 1):
        path = os.path.join(MESH_DIR, fn)
        sid = os.path.splitext(fn)[0]
        grp = "EXP" if sid.startswith("EXP") else "SDG"

        n_raw, fmt = raw_face_count(path)

        m = o3d.io.read_triangle_mesh(path)
        v_raw = len(m.vertices)
        m = m.simplify_quadric_decimation(TARGET_FACES)
        v_dec = np.asarray(m.vertices)
        f_dec = np.asarray(m.triangles)

        f_dec = f_dec[np.all(f_dec < len(v_dec), axis=1)]
        tm = trimesh.Trimesh(vertices=v_dec, faces=f_dec, process=True)
        tm.remove_unreferenced_vertices()
        filter_laplacian(tm, iterations=SMOOTH_ITERS, volume_constraint=False)

        rows.append(dict(
            ID=sid, group=grp, stl_format=fmt,
            faces_raw=n_raw, vertices_raw=v_raw,
            faces_decimated=int(len(tm.faces)),
            vertices_decimated=int(len(tm.vertices)),
            decimation_active=bool(n_raw > TARGET_FACES),
        ))
        if i % 20 == 0 or i == len(files):
            print(f"  {i}/{len(files)}", flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV} ({len(df)} rows)")
    print(df.groupby("group")[["faces_raw", "faces_decimated",
                               "vertices_decimated"]].describe().T.to_string())
    n_inactive = int((~df["decimation_active"]).sum())
    print(f"\nspecimens below the {TARGET_FACES:,}-face target "
          f"(decimation is a no-op): {n_inactive}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
