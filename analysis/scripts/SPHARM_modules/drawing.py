import mesh_processing
import trimesh
import igl
from trimesh.smoothing import filter_laplacian
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from matplotlib.colors import LinearSegmentedColormap


def visualize_mesh(stl_path, target_faces=1000, laplacian_iterations=3):
    """
    Load, simplify, normalize and visualize a mesh.
    Used for visual inspection in SPHARM workflow.
    """
    vertices, faces = mesh_processing.clean_mesh(stl_path)

    v_igl  = np.asfortranarray(np.asarray(vertices), dtype=np.float64)
    f_igl  = np.asfortranarray(np.asarray(faces),    dtype=np.int32)
    result = igl.decimate(v_igl, f_igl, int(target_faces))

    # igl.decimate return (V_out, F_out, face_indices, vertex_indices)
    decimated_vertices = np.array(result[0])
    decimated_faces    = np.array(result[1]).reshape(-1, 3)

    mesh = trimesh.Trimesh(vertices=decimated_vertices,
                           faces=decimated_faces, process=False)
    filter_laplacian(mesh, iterations=laplacian_iterations)
    decimated_vertices = mesh.vertices

    mesh_processing.visualize_error(vertices, decimated_vertices)
    normalized_vertices = mesh_processing.normalize_mesh(decimated_vertices)
    mesh_processing.visualize_with_lighting_and_wire(normalized_vertices,
                                                      decimated_faces)


def visualize_spherical_harmonic_real(l=2, m=2,
                                       resolution_phi=200,
                                       resolution_theta=100,
                                       colors=["#66c2a5", "white", "orange"]):
    """
    Visualize the real part of a spherical harmonic Y_l^m as a 3D petal plot.

    The radial distance at each point is set to |Re(Y_l^m)|,
    and color encodes the sign of Re(Y_l^m).

    Parameters
    ----------
    l : int
        Degree of the spherical harmonic.
    m : int
        Order of the spherical harmonic (|m| <= l).
    resolution_phi : int
        Number of azimuthal sampling points.
    resolution_theta : int
        Number of polar sampling points.
    colors : list of str
        Three colors mapped to [negative, zero, positive] values.
    """
    cmap = LinearSegmentedColormap.from_list(
        "custom", colors, N=256)

    # scipy.special.sph_harm(m, l, azimuth, polar)
    phi,   theta   = np.meshgrid(
        np.linspace(0, 2*np.pi, resolution_phi),
        np.linspace(0,   np.pi, resolution_theta)
    )
    Y_lm  = sph_harm(m, l, phi, theta)
    r     = np.real(Y_lm)
    r_abs = np.abs(r)

    x = r_abs * np.sin(theta) * np.cos(phi)
    y = r_abs * np.sin(theta) * np.sin(phi)
    z = r_abs * np.cos(theta)

    norm         = plt.Normalize(-r.max(), r.max())
    colors_mapped = cmap(norm(r))

    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z,
                    facecolors=colors_mapped,
                    rstride=1, cstride=1,
                    linewidth=0, antialiased=False)
    ax.set_axis_off()
    ax.set_box_aspect([1, 1, 1])

    max_range = r_abs.max()
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])

    plt.title(f"Spherical Harmonic $Y_{{{l}}}^{{{m}}}$ Real Part",
              fontsize=16, pad=20)
    plt.show()


if __name__ == "__main__":
    import os
    base_dir = os.path.dirname(os.path.abspath(__file__))
    stl_file = os.path.abspath(
        os.path.join(base_dir, "..", "3D_models_cores", "55.stl"))

    visualize_mesh(stl_file, target_faces=1000, laplacian_iterations=3)
    visualize_spherical_harmonic_real(l=0, m=0)