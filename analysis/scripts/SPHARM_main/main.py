
import igl
import numpy as np
import pyshtools as pysh
import trimesh
from trimesh.smoothing import filter_laplacian
import os
from SPHARM_modules import mesh_processing, pca_align, reconstruction, spherical_harmonics


def main():

    # 1. 读取文件-read stl files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    stl_path = os.path.join(base_dir, "3d_models_QSY_spheroids_multi_poly", "2_Ellipsoid.stl")
    vertices, faces = mesh_processing.clean_mesh(stl_path)

    # 2. 清理和统一目标面数- clean and optimize mesh
    target_faces = 20000
    result = igl.decimate(vertices, faces, target_faces)
    decimated_vertices = result[1]
    decimated_faces = result[2]
    mesh = trimesh.Trimesh(vertices=decimated_vertices, faces=decimated_faces)
    filter_laplacian(mesh, iterations=3)
    decimated_vertices = mesh.vertices

    # 可视化删减后面片及误差-Visualize optimized mesh and error
    mesh_processing.visualize_error(vertices, decimated_vertices)

    # 3. 执行归一化处理-Normalize mesh to unit sphere
    normalized_vertices = mesh_processing.normalize_mesh(decimated_vertices)

    # 可视化归一化后网格-Visualize normalized mesh
    mesh_processing.visualize_normalization(normalized_vertices, decimated_faces)

    # 方向对齐-alignment
    align_vertices = pca_align.robust_pca_alignment(normalized_vertices, enforce_direction=True, verbose=True)

    # 4. 转换到球坐标-Convert Cartesian coordinates to spherical coordinates
    spherical_coords = spherical_harmonics.cartesian_to_spherical(align_vertices)
    R = spherical_coords[:, 0]
    theta = spherical_coords[:, 1]
    phi = spherical_coords[:, 2]

    # 5. 执行插值-interpolate to regular grid
    grid_size = 256
    grid_r = spherical_harmonics.spherical_interpolate(R, theta, phi, grid_size)
    print(f"r: {np.nanmin(grid_r):.3f} ~ {np.nanmax(grid_r):.3f}")
    spherical_harmonics.visualize_interpolated(grid_r)

    # 6. 计算球谐系数-SPHARM decomposition
    if grid_r is not None:
        clm = spherical_harmonics.compute_spherical_harmonics(
            grid_r,
            normalize=True,
            normalization_method='zero-component'
        )

    # 旋转不变谱-rotation-invariant SPHARM power spectrum
    clm_sh = pysh.SHCoeffs.from_array(clm)
    full, spectrum = spherical_harmonics.process_spherical_harmonics(clm_sh, 'output.csv')
    # print(spectrum.head())
    # print(spectrum.head())
    spherical_harmonics.visualize_power_spectrum(spectrum, max_degree=30, log_scale=True, filename=None)

    # 保存球谐系数-Save spherical harmonic coefficients
    clm_array = spherical_harmonics.clm_to_1d_standard(clm)
    csv_filename = "SHPARM_coeffs_1D.csv"
    np.savetxt(csv_filename, clm_array, delimiter=",", fmt="%.6e")

    # 7. 根据球谐展开重建表面-Surface reconstruction
    clm_sh_truncated = clm_sh.pad(lmax=10)
    reconstructed_grid_truncated = clm_sh_truncated.expand(grid='DH')
    reconstruction.visualize_spherical_harmonics_reconstruction(reconstructed_grid_truncated)


if __name__ == "__main__":
    main()



