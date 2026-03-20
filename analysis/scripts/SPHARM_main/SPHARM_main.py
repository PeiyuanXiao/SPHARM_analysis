import os
import numpy as np
import pandas as pd
import pyshtools as pysh
import igl
import trimesh
from trimesh.smoothing import filter_laplacian
import time
from SPHARM_modules import mesh_processing, pca_align, spherical_harmonics, statistics_analysis
from SPHARM_modules import sphericity as sphericity_module
from SPHARM_modules import curvature as curvature_module

# ============================================================
# Parameter configuration
# ============================================================
TARGET_FACES = 20000
GRID_SIZE    = 256
LMAX         = 20


def process_single_mesh(stl_path):
    """Process a single STL file and return a dictionary containing all features."""
    specimen_id = os.path.splitext(os.path.basename(stl_path))[0]

    # 1. Clean the mesh
    vertices, faces = mesh_processing.clean_mesh(stl_path)
    n_faces = faces.shape[0]
    print(f"  → {specimen_id}: {n_faces:,} faces, decimating to {TARGET_FACES:,}...")

    # 2. Simplify
    v_igl = np.asfortranarray(vertices.astype(np.float64))
    f_igl = np.asfortranarray(faces.astype(np.int32))
    if n_faces <= TARGET_FACES:
        print(f"  → Skipping decimation (faces {n_faces} <= target {TARGET_FACES})")
        decimated_vertices = vertices
        decimated_faces    = faces
    else:
        result = igl.decimate(v_igl, f_igl, int(TARGET_FACES))
        decimated_vertices = np.array(result[0])
        decimated_faces    = np.array(result[1]).reshape(-1, 3)
        if decimated_vertices.shape[0] == 0:
            raise ValueError("igl.decimate returned empty mesh")

    # 3. Smooth
    mesh = trimesh.Trimesh(vertices=decimated_vertices,
                           faces=decimated_faces, process=False)
    filter_laplacian(mesh, iterations=3)
    decimated_vertices = mesh.vertices

    # 4. Normalize + PCA alignment
    normalized_vertices = mesh_processing.normalize_mesh(
                              decimated_vertices, decimated_faces)
    aligned_vertices, _ = pca_align.robust_pca_alignment(
                              normalized_vertices, enforce_direction=True)

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

    # 7. Power spectrum + spectral harmonic entropy (SHE)
    _, spectrum = spherical_harmonics.process_spherical_harmonics(clm_sh)
    total_power = spectrum["total_power"].astype(float)
    SHE         = float(np.sum(clm_sh.spectrum()))

    # 8. Construct result row
    row = {
        "specimen_id":      specimen_id,
        "SHE":              SHE,
        "n_faces_original": int(n_faces),
    }
    for l, p in enumerate(total_power):
        row[f"power_degree_{l}"] = float(p)
    for j, c in enumerate(clm_sh.coeffs.flatten()):
        row[f"coeff_{j}"] = float(np.real(c))

    print(f"  ✓ {specimen_id} (SHE={SHE:.4f})")
    return row


def batch_process(input_dir, output_dir):
    """Batch process all STL files, writing rows to CSV with support for resuming."""
    os.makedirs(output_dir, exist_ok=True)

    output_csv = os.path.join(output_dir, "SPHARM_results.csv")
    failed_csv = os.path.join(output_dir, "SPHARM_failed.csv")

    stl_files = sorted([
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.lower().endswith('.stl')
    ])
    total = len(stl_files)

    # Resume processing from last checkpoint
    processed_ids  = set()
    header_written = False
    if os.path.exists(output_csv):
        df_existing    = pd.read_csv(output_csv)
        processed_ids  = set(df_existing["specimen_id"].astype(str).tolist())
        header_written = True

    failed_ids = set()
    if os.path.exists(failed_csv):
        df_failed  = pd.read_csv(failed_csv)
        failed_ids = set(df_failed["specimen_id"].astype(str).tolist())

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
                failed_rows.append({"specimen_id": specimen_id,
                                    "error": str(e)})
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
    # Post-processing 1: Variance analysis + UMAP
    # --------------------------------------------------------
    print(f"\n{'='*55}")
    print("Running post-processing analysis...")
    print(f"{'='*55}")
    try:
        statistics_analysis.run_batch_analysis(output_csv, output_dir, LMAX)
    except Exception as e:
        print(f"Post-processing failed (non-critical): {e}")

    # --------------------------------------------------------
    # Post-processing 2: Sphericity + centroid offset
    # --------------------------------------------------------
    print(f"\n{'='*55}")
    print("Running sphericity calculation...")
    print(f"{'='*55}")
    sphericity_csv = os.path.join(output_dir, "sphericity_iso.csv")
    try:
        sphericity_module.batch_calculate_sphericity(
            input_dir, sphericity_csv)
    except Exception as e:
        print(f"Sphericity calculation failed (non-critical): {e}")

    # --------------------------------------------------------
    # Post-processing 3: Curvature
    # --------------------------------------------------------
    print(f"\n{'='*55}")
    print("Running curvature calculation...")
    print(f"{'='*55}")
    curvature_csv = os.path.join(output_dir, "curvature.csv")
    try:
        df_curv = curvature_module.batch_average_curvature(input_dir)
        df_curv.to_csv(curvature_csv, index=False)
        print(f"Curvature saved to: {curvature_csv}")
    except Exception as e:
        print(f"Curvature calculation failed (non-critical): {e}")


if __name__ == "__main__":
    input_directory  = "/project/analysis/data/3D_models_cores"
    output_directory = "/project/analysis/data/derived_data"
    batch_process(input_directory, output_directory)
