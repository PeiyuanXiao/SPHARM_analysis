import os
import struct
import tempfile
import gc
import numpy as np
import pandas as pd
import trimesh
import igl
from trimesh.smoothing import filter_laplacian
from SPHARM_modules import mesh_processing

# ============================================================
# [修改] 新增预降采样参数，与 SPHARM_main.py 保持一致
# ============================================================
PRE_DECIMATE_THRESHOLD = 3_000_000
PRE_DECIMATE_TARGET    = 500_000
# ============================================================
# [修改结束]
# ============================================================


def calculate_iso_sphericity(stl_file_path, target_faces=15000):
    """
    Calculate sphericity, surface area and centroid offset for a single STL.

    Uses the simplified and smoothed mesh for all calculations,
    ensuring consistency with the SPHARM main pipeline.

    Parameters
    ----------
    stl_file_path : str
    target_faces : int

    Returns
    -------
    sphericity : float
        ISO sphericity = (π^(1/3) * (6V)^(2/3)) / A
    surface_area_cm2 : float
    centroid_offset : float
        Euclidean distance between mass centroid and bounding box centroid (mm).
    """
    # ============================================================
    # [修改] 预降采样：与 SPHARM_main.py 相同的流式抽样 + open3d 路径，
    #        解决超大文件 OOM 和退化面片导致 igl.decimate 失败的问题。
    # ============================================================
    tmp_path = None
    with open(stl_file_path, 'rb') as f:
        f.read(80)
        n_raw = struct.unpack('<I', f.read(4))[0]

    if n_raw > PRE_DECIMATE_THRESHOLD:
        tmp_path = tempfile.mktemp(suffix='.stl')
        step     = max(1, n_raw // PRE_DECIMATE_TARGET)
        keep_set = set(range(0, n_raw, step))
        n_keep   = len(keep_set)

        with open(stl_file_path, 'rb') as fin, open(tmp_path, 'wb') as fout:
            header = fin.read(80)
            fin.read(4)
            fout.write(header)
            fout.write(struct.pack('<I', n_keep))
            for i in range(n_raw):
                face_data = fin.read(50)
                if i in keep_set:
                    fout.write(face_data)
        gc.collect()
    # ============================================================
    # [修改结束]
    # ============================================================

    try:
        # ============================================================
        # [修改] 超大文件：open3d 路径；普通文件：原有 igl 路径
        # ============================================================
        if tmp_path is not None:
            import open3d as o3d
            o3d_mesh = o3d.io.read_triangle_mesh(tmp_path)
            o3d_mesh = o3d_mesh.simplify_quadric_decimation(target_faces)
            decimated_vertices = np.asarray(o3d_mesh.vertices)
            decimated_faces    = np.asarray(o3d_mesh.triangles)
            if len(decimated_faces) == 0:
                raise ValueError("open3d.simplify_quadric_decimation returned empty mesh")
        else:
            # 1. Clean the mesh
            vertices, faces = mesh_processing.clean_mesh(stl_file_path)

            # 2. Simplify
            v_igl  = np.asfortranarray(vertices.astype(np.float64))
            f_igl  = np.asfortranarray(faces.astype(np.int32))
            result = igl.decimate(v_igl, f_igl, int(target_faces))
            decimated_vertices = np.array(result[0])
            decimated_faces    = np.array(result[1]).reshape(-1, 3)
            if decimated_vertices.shape[0] == 0:
                raise ValueError("igl.decimate returned empty mesh")
        # ============================================================
        # [修改结束]
        # ============================================================

        # 3. Smooth
        valid_mask = np.all(decimated_faces < len(decimated_vertices), axis=1)
        decimated_faces = decimated_faces[valid_mask]

        mesh = trimesh.Trimesh(vertices=decimated_vertices,
                               faces=decimated_faces, process=True)

        if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
            raise ValueError("Mesh became empty after structural cleaning.")

        filter_laplacian(mesh, iterations=3, volume_constraint=False)

        # 4. Watertight repair
        if not mesh.is_watertight:
            mesh.fill_holes()

        # 5. Compute volume and surface area (in mm units)
        volume_mm3       = mesh.volume
        surface_area_mm2 = mesh.area

        if volume_mm3 <= 0:
            raise ValueError(f"Invalid volume: {volume_mm3:.4f} mm³. "
                             f"Check mesh watertightness.")

        # 6. Unit conversion: mm → cm
        volume_cm3       = volume_mm3    / 1000.0
        surface_area_cm2 = surface_area_mm2 / 100.0

        # 7. ISO sphericity formula
        sphericity = ((np.pi ** (1/3)) * (6 * volume_cm3) ** (2/3)) / surface_area_cm2

        # 8. Centroid offset
        centroid_offset = compute_centroid_offset(mesh)

        return sphericity, surface_area_cm2, centroid_offset

    finally:
        # ============================================================
        # [修改] 清理临时文件
        # ============================================================
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)
        # ============================================================
        # [修改结束]
        # ============================================================


def compute_centroid_offset(mesh):
    """
    Compute the Euclidean distance between the mass centroid
    and the bounding box centroid (mm).

    Parameters
    ----------
    mesh : trimesh.Trimesh

    Returns
    -------
    distance : float
    """
    centroid    = mesh.centroid
    bbox_center = mesh.bounding_box.centroid
    distance    = np.linalg.norm(centroid - bbox_center)
    return distance


def batch_calculate_sphericity(input_dir, output_csv, target_faces=15000):
    """
    Batch calculate ISO sphericity for all STL files in a directory.

    Parameters
    ----------
    input_dir : str
    output_csv : str
    target_faces : int
    """
    results   = []
    stl_files = sorted([f for f in os.listdir(input_dir)
                        if f.lower().endswith('.stl')])

    for filename in stl_files:
        stl_path  = os.path.join(input_dir, filename)
        base_name = os.path.splitext(filename)[0]

        # Parse ID and typology from filename
        parts = base_name.split('-')
        if len(parts) >= 2:
            id_part  = "-".join(parts[:-1])
            typology = parts[-1]
        else:
            id_part  = base_name
            typology = ''

        try:
            sphericity, surface_area, centroid_offset = \
                calculate_iso_sphericity(stl_path, target_faces=target_faces)
            results.append({
                'ID':                  id_part,
                'Typology':            typology,
                'sphericity':          round(sphericity,      6),
                'surface_area_cm2':    round(surface_area,    4),
                'centroid_offset_mm':  round(centroid_offset, 6),
                'status':              'success'
            })
            print(f"  ✓ {filename} | "
                  f"sphericity={sphericity:.4f} | "
                  f"area={surface_area:.2f} cm²")

        except Exception as e:
            results.append({
                'ID':                 id_part,
                'Typology':           typology,
                'sphericity':         None,
                'surface_area_cm2':   None,
                'centroid_offset_mm': None,
                'status':             f'error: {str(e)}'
            })
            print(f"  ✗ {filename}: {e}")

    # Sort by ID before saving
    df = pd.DataFrame(results).sort_values('ID').reset_index(drop=True)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)

    success = (df['status'] == 'success').sum()
    print(f"\nFinished! Success: {success}/{len(df)}, "
          f"Failed: {len(df)-success}")

    return df


if __name__ == "__main__":
    base_dir   = os.path.dirname(os.path.abspath(__file__))
    input_dir  = os.path.abspath(
        os.path.join(base_dir, "..", "3D_models_cores"))
    output_csv = os.path.abspath(
        os.path.join(base_dir, "..", "..", "..", "data",
                     "derived_data", "sphericity_iso.csv"))

    batch_calculate_sphericity(input_dir, output_csv)
