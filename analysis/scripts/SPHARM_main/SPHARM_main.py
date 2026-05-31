# Standard library
import os
import sys
import time
import struct
import tempfile
import gc
from pathlib import Path

# Path setup
sys.path.insert(0, str(Path(__file__).parent.parent))

# Third-party
import numpy as np
import pandas as pd
import pyshtools as pysh
import trimesh
from trimesh.smoothing import filter_laplacian

# Local modules
from SPHARM_modules import mesh_processing, pca_align, spherical_harmonics, statistics_analysis
from SPHARM_modules.power_spectrum import compute_power_spectrum

# ============================================================
# Parameter configuration
# ============================================================
TARGET_FACES = 20000
GRID_SIZE    = 256
LMAX         = 20

# ============================================================
# Pre-decimation parameters
#   PRE_DECIMATE_THRESHOLD : trigger pre-decimation when face count exceeds this
#   PRE_DECIMATE_TARGET    : target face count for pre-decimation
#   (after pre-decimation, open3d still reduces to TARGET_FACES)
# ============================================================
PRE_DECIMATE_THRESHOLD = 3_000_000
PRE_DECIMATE_TARGET    = 500_000


def process_single_mesh(stl_path):
    """Process a single STL file and return a dictionary containing all features."""

    import open3d as o3d

    specimen_id = os.path.splitext(os.path.basename(stl_path))[0]

    # ============================================================
    # Use open3d for all files (ASCII and binary) to avoid the
    # methodological inconsistency of mixing igl and open3d.
    # Steps:
    #   1. read the header to detect format and face count
    #   2. above threshold -> pre-decimate to PRE_DECIMATE_TARGET
    #      below threshold -> load directly
    #   3. open3d.simplify_quadric_decimation reduces to TARGET_FACES
    # ============================================================
    tmp_path = None

    with open(stl_path, 'rb') as f:
        header_bytes = f.read(80)
        is_ascii = header_bytes.lstrip().startswith(b'solid')
        if not is_ascii:
            n_raw = struct.unpack('<I', f.read(4))[0]
        else:
            n_raw = None

    if is_ascii:
        print(f"  -> {specimen_id}: ASCII STL detected, loading via open3d...")
        o3d_mesh = o3d.io.read_triangle_mesh(stl_path)
        n_raw = len(o3d_mesh.triangles)
        print(f"  -> {specimen_id}: {n_raw:,} faces")
        if n_raw > PRE_DECIMATE_THRESHOLD:
            print(f"  -> pre-decimating to {PRE_DECIMATE_TARGET:,}...")
            o3d_mesh = o3d_mesh.simplify_quadric_decimation(PRE_DECIMATE_TARGET)
            o3d_mesh.compute_vertex_normals()
            tmp_path = tempfile.mktemp(suffix='.stl')
            o3d.io.write_triangle_mesh(tmp_path, o3d_mesh)
            del o3d_mesh
            gc.collect()
            print(f"  -> pre-decimation done")
            load_path  = tmp_path
            need_final = True
        else:
            load_path     = None
            need_final    = True
            o3d_mesh_hold = o3d_mesh
    else:
        print(f"  -> {specimen_id}: {n_raw:,} faces")
        if n_raw > PRE_DECIMATE_THRESHOLD:
            print(f"  -> pre-decimating to {PRE_DECIMATE_TARGET:,} via streaming...")
            tmp_path = tempfile.mktemp(suffix='.stl')
            step     = max(1, n_raw // PRE_DECIMATE_TARGET)
            keep_set = set(range(0, n_raw, step))
            n_keep   = len(keep_set)
            with open(stl_path, 'rb') as fin, open(tmp_path, 'wb') as fout:
                header = fin.read(80)
                fin.read(4)
                fout.write(header)
                fout.write(struct.pack('<I', n_keep))
                for i in range(n_raw):
                    face_data = fin.read(50)
                    if i in keep_set:
                        fout.write(face_data)
            gc.collect()
            print(f"  -> Streaming done: {n_keep:,} faces saved")
            load_path  = tmp_path
            need_final = True
        else:
            load_path  = stl_path
            need_final = True

    try:
        # All paths reduce to TARGET_FACES here via open3d
        if is_ascii and n_raw <= PRE_DECIMATE_THRESHOLD:
            o3d_mesh_final = o3d_mesh_hold
        else:
            o3d_mesh_final = o3d.io.read_triangle_mesh(load_path)

        o3d_mesh_final = o3d_mesh_final.simplify_quadric_decimation(TARGET_FACES)
        decimated_vertices = np.asarray(o3d_mesh_final.vertices)
        decimated_faces    = np.asarray(o3d_mesh_final.triangles)
        n_faces = len(decimated_faces)
        print(f"  -> {specimen_id}: decimation done: {n_faces:,} faces")

        if n_faces == 0:
            raise ValueError("open3d.simplify_quadric_decimation returned empty mesh")

        # 3. Smooth
        valid_mask = np.all(decimated_faces < len(decimated_vertices), axis=1)
        decimated_faces = decimated_faces[valid_mask]
        mesh = trimesh.Trimesh(vertices=decimated_vertices,
                               faces=decimated_faces, process=True)
        mesh.remove_unreferenced_vertices()
        if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
            raise ValueError("Mesh empty after cleaning unreferenced vertices")
        filter_laplacian(mesh, iterations=3, volume_constraint=False)
        decimated_vertices = mesh.vertices
        decimated_faces    = mesh.faces

        # 4. Normalize + PCA alignment
        normalized_vertices = mesh_processing.normalize_mesh(
                                  decimated_vertices, decimated_faces)
        aligned_vertices, _ = pca_align.robust_pca_alignment(
                                  normalized_vertices,
                                  faces=decimated_faces,
                                  enforce_direction=True)

        # 5. Spherical coordinate interpolation
        spherical_coords = spherical_harmonics.cartesian_to_spherical(aligned_vertices)
        R, theta, phi    = spherical_coords.T
        grid_r           = spherical_harmonics.spherical_interpolate(
                               R, theta, phi, GRID_SIZE)

        # 6. Spherical harmonic expansion
        clm    = spherical_harmonics.compute_spherical_harmonics(
                     grid_r, normalization_method='zero-component')
        clm_sh = pysh.SHCoeffs.from_array(
                     clm, normalization='4pi', csphase=1, lmax=LMAX
                 ).pad(lmax=LMAX)

        # 7. Power spectrum
        norm_power = compute_power_spectrum(clm_sh, lmax=LMAX)["norm_power"]

        # 8. Construct result row
        row = {
            "ID":               specimen_id,
            "n_faces_original": int(n_raw),
        }
        for l, p in enumerate(norm_power):
            row[f"power_l{l}"] = float(p)
        for j, c in enumerate(clm_sh.coeffs.flatten()):
            row[f"coeff_{j}"] = float(np.real(c))

        print(f"  -> {specimen_id}")
        return row

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


def batch_process(input_dir, output_dir):
    """Batch process all STL files, writing rows to CSV with support for resuming."""
    os.makedirs(output_dir, exist_ok=True)

    output_csv = os.path.join(output_dir, "SPHARM_morphology.csv")
    failed_csv = os.path.join(output_dir, "SPHARM_failed.csv")

    stl_files = sorted([
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.lower().endswith('.stl')
    ])
    total = len(stl_files)

    # Resume processing from last checkpoint
    # Use the unified column name "id" when reading processed IDs
    processed_ids  = set()
    header_written = False
    if os.path.exists(output_csv):
        df_existing    = pd.read_csv(output_csv)
        processed_ids  = set(df_existing["ID"].astype(str).tolist())
        header_written = True

    failed_ids = set()
    if os.path.exists(failed_csv):
        df_failed  = pd.read_csv(failed_csv)
        failed_ids = set(df_failed["ID"].astype(str).tolist())

    stl_files_todo = [
        f for f in stl_files
        if os.path.splitext(os.path.basename(f))[0] not in processed_ids
        and os.path.splitext(os.path.basename(f))[0] not in failed_ids
    ]

    print(f"\n{'='*55}")
    print(f"Found {total} STL files total.")
    print(f"  Already processed: {len(processed_ids)}")
    print(f"  Previously failed: {len(failed_ids)}")
    print(f"  To process now:    {len(stl_files_todo)}")
    print(f"{'='*55}\n")

    if not stl_files_todo:
        print("All specimens already processed!")
    else:
        start_time    = time.time()
        success_count = 0
        fail_count    = 0
        failed_rows   = []

        for i, stl_path in enumerate(stl_files_todo, 1):
            specimen_id = os.path.splitext(os.path.basename(stl_path))[0]
            print(f"\n[{i}/{len(stl_files_todo)}] {specimen_id}")

            try:
                row    = process_single_mesh(stl_path)
                df_row = pd.DataFrame([row])
                df_row.to_csv(output_csv, mode='a',
                              header=not header_written, index=False)
                header_written = True
                success_count += 1

            except Exception as e:
                print(f"  ✗ Error: {specimen_id}: {e}")
                import traceback
                traceback.print_exc()
                # failed_csv also uses the unified column name "id"
                failed_rows.append({"ID": specimen_id, "error": str(e)})
                fail_count += 1

        if failed_rows:
            pd.DataFrame(failed_rows).to_csv(failed_csv, index=False)
            print(f"\nFailed specimens saved to: {failed_csv}")

        elapsed = time.time() - start_time
        print(f"\n{'='*55}")
        print(f"SPHARM batch completed!")
        print(f"  Success: {success_count}")
        print(f"  Failed:  {fail_count}")
        print(f"  Time:    {elapsed:.1f}s "
              f"({elapsed/max(success_count, 1):.1f}s per specimen)")
        print(f"  Results: {output_csv}")
        print(f"{'='*55}")

    # --------------------------------------------------------
    # Post-processing 1: Variance analysis
    # --------------------------------------------------------
    print(f"\n{'='*55}")
    print("Running post-processing analysis...")
    print(f"{'='*55}")
    try:
        statistics_analysis.run_batch_analysis(output_csv, output_dir, LMAX)
    except Exception as e:
        print(f"Post-processing failed (non-critical): {e}")

if __name__ == "__main__":
    input_directory  = "/project/analysis/data/3D_models_cores"
    output_directory = "/project/analysis/data/derived_data"
    batch_process(input_directory, output_directory)
