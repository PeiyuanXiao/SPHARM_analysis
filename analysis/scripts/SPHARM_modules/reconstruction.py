import pyvista as pv
import numpy as np


def visualize_spherical_harmonics_reconstruction(grid_sh):
    """
    Visualize a 3D shape reconstructed from spherical harmonics.

    Parameters
    ----------
    grid_sh : pyshtools SHGrid object
        DH grid from pyshtools clm.expand(grid='DH').
        Expected shape: (n_lat, n_lon) where n_lon = 2 * n_lat.
    """
    grid_data = np.real(grid_sh.data)
    n_lat     = grid_data.shape[0]
    n_lon     = grid_data.shape[1]

    theta = np.linspace(0, np.pi,    n_lat, endpoint=True)
    phi   = np.linspace(0, 2*np.pi,  n_lon, endpoint=True)

    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing='ij')
    r = grid_data

    x = (r * np.sin(theta_grid) * np.cos(phi_grid)).T
    y = (r * np.sin(theta_grid) * np.sin(phi_grid)).T
    z = (r * np.cos(theta_grid)).T

    grid = pv.StructuredGrid(x, y, z)

    plotter = pv.Plotter()
    plotter.add_mesh(grid,
                     scalars=r.flatten(),
                     cmap="coolwarm",
                     opacity=1.0,
                     show_edges=False,
                     specular=0.8)
    plotter.add_axes(box_args={'color': 'red'})
    plotter.add_title(f"Reconstruction ({n_lat}×{n_lon})")
    plotter.show()