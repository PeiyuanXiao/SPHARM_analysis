import os
import numpy as np
import pandas as pd
import trimesh
import igl
from trimesh.smoothing import filter_laplacian
from SPHARM_modules import mesh_processing


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
    # 1. Clean the mesh
    vertices, faces = mesh_processing.clean_mesh(stl_file_path)

    # 2. Simplify
    v_igl  = np.asfortranarray(vertices.astype(np.float64))
    f_igl  = np.asfortranarray(faces.astype(np.int32))
    result = igl.decimate(v_igl, f_igl, int(target_faces))
    decimated_vertices = np.array(result[0])
    decimated_faces    = np.array(result[1]).reshape(-1, 3)

    # 3. Smooth
    mesh = trimesh.Trimesh(vertices=decimated_vertices,
                           faces=decimated_faces, process=False)
    filter_laplacian(mesh, iterations=3)

    # 4. Watertight repair
    if not mesh.is_watertight:
        mesh.fill_holes()

    # 5. Compute volume and surface area (in mm units)
    volume_mm3     = mesh.volume
    surface_area_mm2 = mesh.area

    if volume_mm3 <= 0:
        raise ValueError(f"Invalid volume: {volume_mm3:.4f} mm³. "
                         f"Check mesh watertightness.")

    # 6. Unit conversion: mm → cm
    volume_cm3       = volume_mm3   / 1000.0
    surface_area_cm2 = surface_area_mm2 / 100.0

    # 7. ISO sphericity formula
    sphericity = ((np.pi ** (1/3)) * (6 * volume_cm3) ** (2/3)) / surface_area_cm2

    # 8. Centroid offset
    centroid_offset = compute_centroid_offset(mesh)

    return sphericity, surface_area_cm2, centroid_offset


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
    results      = []
    stl_files    = sorted([f for f in os.listdir(input_dir)
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
                calculate_iso_sphericity(stl_path,
                                         target_faces=target_faces)
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