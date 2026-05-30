# pca_align.py
# Area-weighted PCA alignment of a 3D point cloud.
# robust_pca_alignment() is called by SPHARM_main.py.

import numpy as np


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
