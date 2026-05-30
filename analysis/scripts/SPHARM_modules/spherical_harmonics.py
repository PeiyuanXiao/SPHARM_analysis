# spherical_harmonics.py
# Cartesian->spherical conversion, grid interpolation, and SH expansion.
# Called by SPHARM_main.py.

import numpy as np
import pyshtools as pysh
import pyshtools.expand as shtools


def cartesian_to_spherical(normalized_vertices):
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    normalized_vertices : np.ndarray, shape (N, 3)

    Returns
    -------
    spherical_coords : np.ndarray, shape (N, 3)
        Columns: [r, theta, phi]
        r     : radial distance
        theta : colatitude in [0, π]
        phi   : longitude  in [0, 2π)
    """
    x, y, z = normalized_vertices.T
    r       = np.sqrt(x**2 + y**2 + z**2)
    r       = np.where(r == 0, 1e-5, r)
    theta   = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi     = np.arctan2(y, x) % (2 * np.pi)
    return np.column_stack([r, theta, phi])


def spherical_interpolate(R, theta, phi, grid_size):
    """
    Interpolate scattered spherical data onto a regular grid.

    Parameters
    ----------
    R : array-like
        Radii values.
    theta : array-like
        Colatitude angles in [0, π].
    phi : array-like
        Longitude angles in [0, 2π).
    grid_size : int
        Number of latitude points. Output shape: (grid_size, grid_size).

    Returns
    -------
    grid : np.ndarray, shape (grid_size, grid_size) or None
    """
    from scipy.interpolate import griddata

    if len(R) < 4:
        return None

    I = np.linspace(0, np.pi,    grid_size, endpoint=False)
    J = np.linspace(0, 2*np.pi,  grid_size, endpoint=False)
    J, I = np.meshgrid(J, I)

    values = R
    points = np.array([theta, phi]).T

    # Polar completion
    points = np.concatenate((points,
                              np.array([[0, 0],
                                        [0, 2*np.pi],
                                        [np.pi, 0],
                                        [np.pi, 2*np.pi]])), axis=0)
    rmin   = np.mean(R[theta == theta.min()])
    rmax   = np.mean(R[theta == theta.max()])
    values = np.concatenate((values, [rmin, rmin, rmax, rmax]))

    # Longitude periodicity handling
    points = np.concatenate((points,
                              points - [0, 2*np.pi],
                              points + [0, 2*np.pi]), axis=0)
    values = np.concatenate((values, values, values))

    xi   = np.array([[I[i, j], J[i, j]]
                     for i in range(grid_size)
                     for j in range(grid_size)])
    grid = griddata(points, values, xi, method='linear')
    grid = grid.reshape((grid_size, grid_size))
    grid[:, -1] = grid[:, 0]

    return grid


def compute_spherical_harmonics(surface,
                                 normalize=True,
                                 normalization_method='zero-component'):
    """
    Compute spherical harmonic coefficients from a 2D surface grid
    using the Driscoll-Healy sampling theorem.

    Parameters
    ----------
    surface : np.ndarray, shape (N, N) or (N, 2N)
        Both dimensions must be even.
    normalize : bool
        Whether to normalize the coefficients.
    normalization_method : {'zero-component', 'mean-radius'}
        - 'zero-component': divide all coefficients by c(l=0, m=0),
          making coefficients dimensionless relative to mean radius.
        - 'mean-radius': divide surface by its mean absolute value
          before expansion.

    Returns
    -------
    harmonics : np.ndarray, shape (2, lmax+1, lmax+1)
        Complex coefficients in pyshtools '4pi' normalization convention,
        as returned by SHExpandDHC. Use with:
            pysh.SHCoeffs.from_array(harmonics, normalization='4pi')

    Notes
    -----
    SHExpandDHC always returns coefficients in '4pi' normalization.
    The zero-component normalization scales the array but does not
    change the underlying convention.
    """
    if surface.shape[1] % 2 or surface.shape[0] % 2:
        raise ValueError("Grid dimensions must be even")

    if surface.shape[1] == surface.shape[0]:
        sampling = 1
    elif surface.shape[1] == 2 * surface.shape[0]:
        sampling = 2
    else:
        raise ValueError("Grid must be (N, N) or (N, 2N)")

    processed_surface = surface.copy()
    if normalize and normalization_method == 'mean-radius':
        processed_surface /= np.mean(np.abs(processed_surface))

    # SHExpandDHC returns coefficients under the '4pi' normalization convention
    harmonics = shtools.SHExpandDHC(processed_surface, sampling=sampling)

    if normalize and normalization_method == 'zero-component':
        c00 = harmonics[0, 0, 0]
        if np.abs(c00) < 1e-10:
            raise ValueError("c(l=0, m=0) is near zero, cannot normalize")
        harmonics = harmonics / c00

    return harmonics
