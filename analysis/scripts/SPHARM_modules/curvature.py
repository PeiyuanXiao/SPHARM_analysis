import trimesh
import igl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os


def simplify_mesh(vertices, faces, target_faces):
    """
    Robustly simplify the mesh using libigl.
    """
    if target_faces >= len(faces):
        return vertices, faces

    try:
        v_igl  = np.asfortranarray(np.asarray(vertices), dtype=np.float64)
        f_igl  = np.asfortranarray(np.asarray(faces),    dtype=np.int32)
        target = int(target_faces)

        res = igl.decimate(v_igl, f_igl, target)

        # igl.decimate return (V_out, F_out, face_indices, vertex_indices)
        if isinstance(res[0], bool) or isinstance(res[0], np.bool_):
            success, simplified_vertices, simplified_faces = res[0], res[1], res[2]
            if not success:
                print("Warning: igl.decimate indicated failure. Returning original.")
                return vertices, faces
        else:
            simplified_vertices, simplified_faces = res[0], res[1]

        if simplified_vertices is None or len(simplified_faces) == 0:
            return vertices, faces

        return simplified_vertices, simplified_faces

    except Exception as e:
        print(f"Decimation error: {e}. Falling back to original mesh.")
        return vertices, faces


def compute_curvature(mesh, k_neighbors=30):
    """
    Compute an approximate curvature for each triangle in a mesh.

    The curvature is estimated as the angle between the triangle's normal
    and the normal of a locally fitted plane through its neighboring triangle
    centers (SVD-based).

    Parameters
    ----------
    mesh : trimesh.Trimesh
    k_neighbors : int
        Number of nearest neighbor triangles used to fit the local plane.
        Note: this value is scale-dependent — the same k_neighbors corresponds
        to different physical areas at different mesh densities. For meaningful
        cross-specimen comparison, ensure all meshes are simplified to the same
        target face count before calling this function.

    Returns
    -------
    curvatures : np.ndarray
        Array of curvature angles (degrees) for each triangle.
    """
    tri_centers = mesh.triangles_center
    tri_normals = mesh.face_normals

    nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1).fit(tri_centers)
    _, indices = nbrs.kneighbors(tri_centers)

    curvatures = []
    for i, neighbors in enumerate(indices):
        if len(neighbors) <= 1:
            curvatures.append(np.nan)
            continue

        neighbor_points = tri_centers[neighbors[1:]]
        _, _, Vt        = np.linalg.svd(
                              neighbor_points - neighbor_points.mean(axis=0))
        local_normal    = Vt[-1]

        dot   = np.clip(np.dot(tri_normals[i], local_normal), -1.0, 1.0)
        angle = np.degrees(np.arccos(abs(dot)))
        curvatures.append(angle)

    return np.array(curvatures)


def batch_average_curvature(input_folder, target_faces=10000, k_neighbors=30):
    """
    Batch compute average curvature for all STL files in a folder.

    Parameters
    ----------
    input_folder : str
    target_faces : int
        All meshes are simplified to this face count before curvature
        computation, ensuring k_neighbors is spatially comparable across
        specimens.
    k_neighbors : int
        See compute_curvature.
    """
    files   = [f for f in os.listdir(input_folder)
               if f.lower().endswith('.stl')]
    results = []

    for file in files:
        path        = os.path.join(input_folder, file)
        mesh        = trimesh.load(path)
        simple_v, simple_f = simplify_mesh(mesh.vertices, mesh.faces,
                                           target_faces)
        simple_mesh = trimesh.Trimesh(vertices=simple_v, faces=simple_f)
        curvatures  = compute_curvature(simple_mesh, k_neighbors)
        avg_curv    = np.nanmean(curvatures)

        filename_no_ext = os.path.splitext(file)[0]
        results.append({
            'filename':             filename_no_ext,
            'average_curvature_deg': avg_curv
        })
        print(f"{filename_no_ext}: average curvature = {avg_curv:.3f} degrees")

    return pd.DataFrame(results)


def visualize_curvature(mesh, curvatures, title="Curvature (degrees)"):
    """Visualize curvature distribution on a mesh."""
    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection='3d')

    norm   = plt.Normalize(vmin=np.nanmin(curvatures),
                           vmax=np.nanmax(curvatures))
    colors = plt.cm.viridis(norm(curvatures))

    collection = Poly3DCollection(mesh.vertices[mesh.faces],
                                  facecolors=colors, edgecolor='none')
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(mesh.vertices[:, 0],
                      mesh.vertices[:, 1],
                      mesh.vertices[:, 2])
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    plt.title(title)

    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    mappable.set_array(curvatures)
    fig.colorbar(mappable, ax=ax, shrink=0.6, aspect=10,
                 label='Curvature (°)')
    plt.show()


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # visualization of single specimen
    stl_file = os.path.abspath(
        os.path.join(base_dir, "..", "3D_models_cores", "55.stl"))

    mesh = trimesh.load(stl_file)
    simple_v, simple_f = simplify_mesh(mesh.vertices, mesh.faces,
                                       target_faces=10000)
    simple_mesh = trimesh.Trimesh(vertices=simple_v, faces=simple_f)
    curvatures  = compute_curvature(simple_mesh)
    visualize_curvature(simple_mesh, curvatures, title="55_Curvature")

    # batch processing
    input_folder  = os.path.abspath(
        os.path.join(base_dir, "..", "3D_models_cores"))
    output_folder = os.path.abspath(
        os.path.join(base_dir, "..", "..", "..", "data", "drived_data"))
    os.makedirs(output_folder, exist_ok=True)

    df = batch_average_curvature(input_folder)
    df.to_csv(os.path.join(output_folder, "curvature.csv"), index=False)
    print(f"Curvature saved to {output_folder}/curvature.csv")