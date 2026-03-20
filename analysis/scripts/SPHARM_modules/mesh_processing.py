import pyvista as pv
import numpy as np
import trimesh
import igl
import pymeshfix
from scipy.spatial import KDTree


def clean_mesh(filepath):
    """
    Clean a 3D mesh from an STL file.

    Steps:
    - Remove unreferenced vertices
    - Fix face index offset if needed
    - Remove infinite values
    - Remove degenerate faces
    - Fix normal orientation and winding consistency

    Parameters
    ----------
    filepath : str

    Returns
    -------
    vertices : np.ndarray
    faces : np.ndarray
    """
    v, f = igl.read_triangle_mesh(filepath)
    v, f, _, _ = igl.remove_unreferenced(v, f)

    if f.min() == 1:
        f = f - 1

    mesh = trimesh.Trimesh(vertices=v, faces=f)
    mesh.remove_infinite_values()

    # Remove degenerate faces
    mask = mesh.nondegenerate_faces().astype(bool)
    mesh.update_faces(mask)

    # Fix normal consistency
    trimesh.repair.fix_normals(mesh)
    trimesh.repair.fix_winding(mesh)
    trimesh.repair.fill_holes(mesh)

    return mesh.vertices, mesh.faces


def normalize_mesh(vertices, faces=None):
    """
    Normalize a mesh so that it is centered at the volume centroid
    and scaled to fit inside a unit sphere.

    Parameters
    ----------
    vertices : np.ndarray, shape (N, 3)
    faces : np.ndarray, shape (M, 3), optional
        If provided, uses volume centroid (mass-weighted) for centering.
        If None, falls back to vertex mean.

    Returns
    -------
    normalized_vertices : np.ndarray, shape (N, 3)
    """
    if faces is not None:
        # Volume centroid
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        centroid = mesh.center_mass
    else:
        # Fallback to vertex mean
        centroid = np.mean(vertices, axis=0)

    centered_vertices = vertices - centroid
    max_radius        = np.max(np.linalg.norm(centered_vertices, axis=1))
    normalized_vertices = centered_vertices / max_radius

    return normalized_vertices


def hausdorff_distance(points1, points2):
    """Compute the Hausdorff distance between two point clouds."""
    d1 = KDTree(points1).query(points2)[0]
    d2 = KDTree(points2).query(points1)[0]
    return max(np.max(d1), np.max(d2))


def visualize_error(vertices, decimated_vertices):
    """Visualize the Hausdorff error between original and decimated mesh."""
    error_mesh = pv.PolyData(vertices)
    error_mesh["Distance"] = KDTree(decimated_vertices).query(vertices)[0]

    plotter = pv.Plotter()
    plotter.add_mesh(error_mesh,
                     scalars="Distance",
                     cmap="coolwarm",
                     opacity=1.0,
                     show_edges=True,
                     scalar_bar_args={"title": "error(mm)"})
    plotter.add_mesh(pv.PolyData(decimated_vertices),
                     color="cyan", show_edges=True, opacity=0.8)
    plotter.show()


def visualize_normalization(normalized_vertices, decimated_faces):
    """Visualize and validate the normalized mesh."""
    pyvista_faces = np.insert(
        decimated_faces.astype(np.int64), 0, 3, axis=1).ravel()

    print("Face array validation:")
    print(f"  Original shape:  {decimated_faces.shape} (expected [n_faces, 3])")
    print(f"  Converted shape: {pyvista_faces.shape} (expected [n_faces*4,])")
    print(f"  Index range:     {decimated_faces.min()} ~ {decimated_faces.max()} "
          f"(should be < {len(normalized_vertices)})")

    normalized_mesh = pv.PolyData(normalized_vertices, pyvista_faces)
    plotter = pv.Plotter()
    plotter.add_mesh(normalized_mesh,
                     color="lightblue",
                     show_edges=True,
                     edge_color="gray",
                     opacity=1,
                     label=f"Mesh after normalization "
                           f"({len(pyvista_faces) // 4} faces)")
    plotter.add_axes(box_args={'color': 'red'})
    plotter.add_title(f"\nvertices: {len(normalized_vertices)}"
                      f"\nfaces: {len(pyvista_faces) // 4}")
    plotter.show()


def visualize_with_lighting_and_wire(normalized_vertices, decimated_faces):
    """Visualization used in SPHARM cluster figures."""
    pyvista_faces = np.insert(
        decimated_faces.astype(np.int64), 0, 3, axis=1).ravel()
    mesh = pv.PolyData(normalized_vertices, pyvista_faces)

    plotter = pv.Plotter(window_size=[800, 800])
    plotter.set_background('white')
    plotter.add_mesh(mesh,
                     color='#f2f2f2',
                     opacity=0.5,
                     show_edges=False,
                     lighting=True,
                     smooth_shading=True)
    plotter.add_mesh(mesh,
                     style='wireframe',
                     color='#999999',
                     line_width=1.2,
                     lighting=False)
    plotter.add_axes()
    plotter.show()