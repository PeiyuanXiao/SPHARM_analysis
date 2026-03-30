import trimesh
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import open3d as o3d
import igl
import copy
import os
from SPHARM_modules import mesh_processing


def clean_mesh(filepath):
    """Clean a 3D mesh from an STL file."""
    v, f = igl.read_triangle_mesh(filepath)
    v, f, _, _ = igl.remove_unreferenced(v, f)

    if f.min() == 1:
        f = f - 1

    mesh = trimesh.Trimesh(vertices=v, faces=f)
    mesh.remove_infinite_values()

    non_deg_mask = mesh.nondegenerate_faces()
    mesh.update_faces(non_deg_mask)

    trimesh.repair.fix_normals(mesh)
    trimesh.repair.fix_winding(mesh)

    return mesh.vertices, mesh.faces


def decimate_mesh(vertices, faces, target_faces=20000):
    """Simplify mesh using libigl."""
    v_igl  = np.asfortranarray(np.asarray(vertices), dtype=np.float64)
    f_igl  = np.asfortranarray(np.asarray(faces),    dtype=np.int32)
    result = igl.decimate(v_igl, f_igl, int(target_faces))

    # igl.decimate return (V_out, F_out, face_indices, vertex_indices)
    v_decim = np.array(result[0])
    f_decim = np.array(result[1]).reshape(-1, 3)

    return v_decim, f_decim


def normalize_mesh(vertices, faces=None):
    """
    Normalize mesh to unit sphere centered at volume centroid.

    Parameters
    ----------
    vertices : np.ndarray, shape (N, 3)
    faces : np.ndarray, shape (M, 3), optional
        If provided, uses volume centroid. If None, falls back to vertex mean.
    """
    if faces is not None:
        mesh     = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        centroid = mesh.center_mass
    else:
        centroid = np.mean(vertices, axis=0)

    centered   = vertices - centroid
    max_radius = np.max(np.linalg.norm(centered, axis=1))
    return centered / max_radius


def robust_pca_alignment(points, faces=None, enforce_direction=True, verbose=False):
    """
    PCA alignment for 3D point clouds, with optional face-area weighting.

    Aligns the input point cloud so that the three principal axes
    correspond to X, Y, Z respectively.

    Sign convention (enforce_direction=True):
        For each axis, the positive direction is defined as the side whose
        points carry more squared projection energy (sum of proj²). This is
        more stable than a median-based rule for near-symmetric shapes,
        because it focuses on the extremes of the distribution rather than
        its centre.

    Area weighting (faces is not None):
        Each vertex is weighted by one-third of the total area of its
        adjacent faces, so high-curvature regions (which are over-sampled
        after decimation) do not bias the principal axes.

    Parameters
    ----------
    points : np.ndarray, shape (N, 3)
    faces : np.ndarray, shape (M, 3), optional
        Triangle indices. If provided, uses face-area weights for PCA.
        If None, falls back to uniform vertex weights.
    enforce_direction : bool
        If True, enforce consistent sign for all three principal axes.
    verbose : bool
        If True, print alignment diagnostics.

    Returns
    -------
    aligned_points : np.ndarray, shape (N, 3)
    rotation_matrix : np.ndarray, shape (3, 3)
    """
    if not isinstance(points, np.ndarray) or points.shape[1] != 3:
        raise ValueError("Input must be an Nx3 NumPy array")
    if len(points) < 3:
        raise ValueError("At least 3 points are required")

    # ── Area-weighted centroid and covariance ──────────────────────────
    if faces is not None:
        # Compute per-face area and distribute 1/3 to each vertex
        v0 = points[faces[:, 0]]
        v1 = points[faces[:, 1]]
        v2 = points[faces[:, 2]]
        face_areas   = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
        vertex_weights = np.zeros(len(points))
        for k in range(3):
            np.add.at(vertex_weights, faces[:, k], face_areas / 3.0)
        total = vertex_weights.sum()
        if total < 1e-12:
            vertex_weights = np.ones(len(points))   # fallback
        else:
            vertex_weights /= total
    else:
        vertex_weights = np.ones(len(points)) / len(points)

    # Weighted centroid
    centroid = (points * vertex_weights[:, None]).sum(axis=0)
    centered = points - centroid

    # Weighted covariance matrix
    cov_matrix = (centered * vertex_weights[:, None]).T @ centered

    _, _, Vt        = np.linalg.svd(cov_matrix)
    rotation_matrix = Vt.T                          # columns = principal axes

    # Ensure right-handed coordinate system
    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix[:, 2] *= -1

    # ── Stable sign convention: energy-based ──────────────────────────
    # For each axis the positive direction is the side that carries more
    # squared-projection energy (weighted). This is robust for near-
    # symmetric shapes where median-based rules can flip arbitrarily.
    if enforce_direction:
        projected = centered @ rotation_matrix
        for axis in range(3):
            proj   = projected[:, axis]
            w      = vertex_weights
            pos_e  = np.dot(w[proj > 0], proj[proj > 0] ** 2) if np.any(proj > 0) else 0.0
            neg_e  = np.dot(w[proj < 0], proj[proj < 0] ** 2) if np.any(proj < 0) else 0.0
            if neg_e > pos_e:
                rotation_matrix[:, axis] *= -1
        # Re-enforce right-handedness after potential flips
        if np.linalg.det(rotation_matrix) < 0:
            rotation_matrix[:, 2] *= -1

    aligned_points = centered @ rotation_matrix

    if verbose:
        print("\n===== PCA Alignment Check =====")
        print(f"Area weighting: {'yes' if faces is not None else 'no (uniform)'}")
        print(f"Determinant of rotation matrix: {np.linalg.det(rotation_matrix):.6f}")
        print(f"PC1 (X): {rotation_matrix[:, 0]}")
        print(f"PC2 (Y): {rotation_matrix[:, 1]}")
        print(f"PC3 (Z): {rotation_matrix[:, 2]}")
        for i, axis in enumerate(['X', 'Y', 'Z']):
            print(f"{axis} range: [{aligned_points[:, i].min():.3f}, "
                  f"{aligned_points[:, i].max():.3f}]")

    return aligned_points, rotation_matrix


def plot_pca_aligned_points(points):
    """Visualize the aligned point cloud with XYZ axes."""
    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection='3d')

    ax.scatter(points[:, 0], points[:, 1], points[:, 2],
               s=1, c='gray', alpha=0.6, label='Point Cloud')

    centroid    = np.mean(points, axis=0)
    max_range   = np.max(np.ptp(points, axis=0))
    axis_length = max_range * 1.1

    for vec, color, label in zip(
        [[1,0,0],[0,1,0],[0,0,1]],
        ['r','g','b'],
        ['X axis','Y axis','Z axis']
    ):
        ax.plot(
            [centroid[0] - axis_length/2 * vec[0],
             centroid[0] + axis_length/2 * vec[0]],
            [centroid[1] - axis_length/2 * vec[1],
             centroid[1] + axis_length/2 * vec[1]],
            [centroid[2] - axis_length/2 * vec[2],
             centroid[2] + axis_length/2 * vec[2]],
            color=color, lw=2, label=label
        )

    ax.scatter(*centroid, s=50, c='k', marker='o', label='Centroid')
    ax.set_box_aspect([1, 1, 1])
    ax.set_title("PCA Aligned Point Cloud")
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.legend()
    plt.show()


def load_template(template_path):
    """Prepare a template model for ICP registration."""
    v_temp, f_temp       = clean_mesh(template_path)
    v_decim, f_decim     = decimate_mesh(v_temp, f_temp, 20000)
    template_normalized  = normalize_mesh(v_decim, f_decim)
    template_aligned, _  = robust_pca_alignment(template_normalized)
    return template_aligned


def prepare_pointcloud(points):
    """Convert a NumPy array to an Open3D PointCloud object."""
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points)
    return pcd


def multi_resolution_icp(source, target,
                          init_transform=np.eye(4), verbose=True):
    """
    Multi-resolution ICP registration.

    Parameters
    ----------
    source : o3d.geometry.PointCloud
    target : o3d.geometry.PointCloud
    init_transform : np.ndarray, shape (4, 4)
    verbose : bool

    Returns
    -------
    transformation : np.ndarray, shape (4, 4)
    """
    voxel_sizes    = [0.1, 0.05, 0.02]
    max_iterations = [200, 100, 50]

    current_transform = init_transform
    final_result      = None

    for voxel_size, max_iter in zip(voxel_sizes, max_iterations):
        source_down = source.voxel_down_sample(voxel_size)
        target_down = target.voxel_down_sample(voxel_size)
        source_down.estimate_normals()
        target_down.estimate_normals()

        result = o3d.pipelines.registration.registration_icp(
            source_down, target_down,
            max_correspondence_distance=voxel_size * 2,
            init=current_transform,
            estimation_method=o3d.pipelines.registration
                               .TransformationEstimationPointToPlane(),
            criteria=o3d.pipelines.registration.ICPConvergenceCriteria(
                max_iteration=max_iter))

        current_transform = result.transformation

        final_result = o3d.pipelines.registration.registration_icp(
            source, target,
            max_correspondence_distance=voxel_sizes[-1],
            init=current_transform,
            estimation_method=o3d.pipelines.registration
                               .TransformationEstimationPointToPlane()
        )

        if verbose:
            evaluate_icp_result(final_result)

    return final_result.transformation


def evaluate_icp_result(result):
    """Print ICP registration diagnostics."""
    try:
        print(f"Fitness:     {result.fitness:.4f} (1 is best)")
        print(f"Inlier RMSE: {result.inlier_rmse * 1000:.2f} mm")

        T    = result.transformation
        R, t = T[:3, :3], T[:3, 3]
        det  = np.linalg.det(R)

        if abs(det - 1) > 1e-3:
            print(f"Warning: rotation matrix determinant = {det:.6f}")

        trans_norm = np.linalg.norm(t)
        print(f"Translation magnitude: {trans_norm * 1000:.2f} mm")
        print(f"Translation vector:    "
              f"[{t[0]*1000:.2f}, {t[1]*1000:.2f}, {t[2]*1000:.2f}] mm")

        cond_num = np.linalg.cond(T)
        print(f"Condition number: {cond_num:.2e}")
        if cond_num > 1e6:
            print("Warning: near-singular matrix, result may be unreliable")

    except Exception as e:
        print(f"Evaluation error: {e}")


def align_all_to_template(template_path, stone_paths):
    """Alignment workflow for multiple models to a template."""
    template_points = load_template(template_path)
    template_pcd    = prepare_pointcloud(template_points)
    template_pcd.estimate_normals()

    aligned_results = []

    for path in stone_paths:
        v_clean, f_clean   = clean_mesh(path)
        v_decim, f_decim   = decimate_mesh(v_clean, f_clean, 20000)
        points_norm        = normalize_mesh(v_decim, f_decim)
        points_pca, _      = robust_pca_alignment(points_norm)

        stone_pcd = prepare_pointcloud(points_pca)
        stone_pcd.estimate_normals()

        transform = multi_resolution_icp(stone_pcd, template_pcd)
        visualize_alignment(stone_pcd, template_pcd, transform)
        stone_pcd.transform(transform)

        aligned_results.append(np.asarray(stone_pcd.points))

    return aligned_results


def visualize_alignment(source_pcd, target_pcd, transform):
    """Visualize before/after alignment."""
    source_original = copy.deepcopy(source_pcd)
    source_aligned  = copy.deepcopy(source_pcd)
    source_aligned.transform(transform)

    source_original.paint_uniform_color([1, 0, 0])
    source_aligned.paint_uniform_color( [0, 0, 1])
    target_pcd.paint_uniform_color(     [0, 1, 0])

    o3d.visualization.draw_geometries(
        [source_original, source_aligned, target_pcd],
        window_name="Alignment Comparison",
        width=1200, height=800, left=200, top=200
    )


def visualize_icp_result(source, target, transform):
    """Visualize ICP result."""
    source_temp = copy.deepcopy(source)
    target_temp = copy.deepcopy(target)
    source_temp.transform(transform)

    source_temp.paint_uniform_color([1, 0.706, 0])
    target_temp.paint_uniform_color([0, 0.651, 0.929])

    coord_frame = o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1)
    o3d.visualization.draw_geometries(
        [source_temp, target_temp, coord_frame],
        window_name="ICP Result",
        width=1200, height=800,
        point_show_normal=True
    )


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    stl_file = os.path.abspath(
        os.path.join(base_dir, "..", "3D_models_cores", "55.stl"))

    vertices, faces   = clean_mesh(stl_file)
    v_decim, f_decim  = decimate_mesh(vertices, faces, target_faces=20000)
    points_norm       = normalize_mesh(v_decim, f_decim)
    aligned, rot_mat  = robust_pca_alignment(points_norm,
                                             faces=f_decim,
                                             verbose=True)
    plot_pca_aligned_points(aligned)
