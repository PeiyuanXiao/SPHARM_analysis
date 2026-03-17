import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
import pyshtools.expand as shtools
import pyshtools as pysh
from scipy.stats import entropy as scipy_entropy


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


def visualize_interpolated(grid_r):
    """Visualize the interpolated regular spherical grid."""
    grid_size = grid_r.shape[0]

    theta = np.linspace(0, np.pi,   grid_size, endpoint=True)
    phi   = np.linspace(0, 2*np.pi, grid_size, endpoint=True)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing='ij')

    print("Sampling check:")
    print(f"  theta: {theta[0]:.3f} ~ {theta[-1]:.3f} (expected 0 ~ π)")
    print(f"  phi:   {phi[0]:.3f} ~ {phi[-1]:.3f} (expected 0 ~ 2π)")

    x = grid_r * np.sin(theta_grid) * np.cos(phi_grid)
    y = grid_r * np.sin(theta_grid) * np.sin(phi_grid)
    z = grid_r * np.cos(theta_grid)

    # Longitude closure check
    tol   = 1e-3
    x_diff = np.max(np.abs(x[:, 0] - x[:, -1]))
    y_diff = np.max(np.abs(y[:, 0] - y[:, -1]))
    z_diff = np.max(np.abs(z[:, 0] - z[:, -1]))
    print("Closure check (φ=0 vs φ=2π):")
    print(f"  x={x_diff:.3e}, y={y_diff:.3e}, z={z_diff:.3e}")
    if max(x_diff, y_diff, z_diff) > tol:
        print("  Warning: coordinates at φ=0 and φ=2π do not match!")

    grid_pv         = pv.StructuredGrid(x, y, z)
    grid_pv["Radius"] = grid_r.T.ravel()
    plotter = pv.Plotter()
    plotter.add_mesh(grid_pv,
                     color='cyan',
                     opacity=0.8,
                     show_edges=True,
                     scalars=grid_r.flatten(),
                     cmap='viridis')
    plotter.add_title(f"Interpolated Grid ({grid_size}×{grid_size})")
    plotter.show()


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


def clm_to_1d_standard(clm, target_l_max=30):
    """
    Convert spherical harmonic coefficients to a 1D real-valued array.

    Parameters
    ----------
    clm : np.ndarray, shape (2, lmax+1, lmax+1)
    target_l_max : int

    Returns
    -------
    np.ndarray, shape ((target_l_max+1)**2,)
        Ordered as: for each l: [c(l,-l), ..., c(l,0), ..., c(l,l)]
    """
    if clm.shape[1] - 1 < target_l_max:
        raise ValueError("Input lmax is smaller than target_l_max")

    coeffs = []
    for l in range(target_l_max + 1):
        for m in range(l, 0, -1):
            coeffs.append(clm[1, l, m])
        coeffs.append(clm[0, l, 0])
        for m in range(1, l + 1):
            coeffs.append(clm[0, l, m])

    return np.real(np.array(coeffs))


def process_spherical_harmonics(clm, output_path=None):
    """
    Analyze spherical harmonic coefficients and compute the power spectrum.

    Parameters
    ----------
    clm : pysh.SHCoeffs or np.ndarray, shape (2, lmax+1, lmax+1)
    output_path : str or None
        If provided, saves results to .xlsx or .csv.

    Returns
    -------
    full_df : pd.DataFrame
        Complete coefficients with degree, order, amplitude, power,
        real, imaginary parts.
    spectrum_df : pd.DataFrame
        Power spectrum aggregated by degree.
    """
    if isinstance(clm, pysh.SHCoeffs):
        coeffs = clm.to_array()
        lmax   = clm.lmax
    elif isinstance(clm, np.ndarray):
        if clm.ndim != 3 or clm.shape[0] != 2:
            raise ValueError("Array must have shape (2, lmax+1, lmax+1)")
        coeffs = clm
        lmax   = clm.shape[1] - 1
    else:
        raise TypeError("Input must be SHCoeffs or ndarray")

    # Real coefficients are automatically converted, compatible with SHExpandDH
    if not np.iscomplexobj(coeffs):
        coeffs = coeffs.astype(complex)

    n_records = (lmax + 1) ** 2
    data = {
        'degree': np.empty(n_records, dtype=np.int32),
        'order':  np.empty(n_records, dtype=np.int32),
        'value':  np.empty(n_records, dtype=np.complex128)
    }

    idx = 0
    for l in range(lmax + 1):
        if l == 0:
            data['degree'][idx] = 0
            data['order'][idx]  = 0
            data['value'][idx]  = coeffs[0, 0, 0]
            idx += 1
        else:
            for m in range(l, 0, -1):
                data['degree'][idx] = l
                data['order'][idx]  = -m
                data['value'][idx]  = coeffs[1, l, m]
                idx += 1
            data['degree'][idx] = l
            data['order'][idx]  = 0
            data['value'][idx]  = coeffs[0, l, 0]
            idx += 1
            for m in range(1, l + 1):
                data['degree'][idx] = l
                data['order'][idx]  = m
                data['value'][idx]  = coeffs[0, l, m]
                idx += 1

    full_df              = pd.DataFrame(data)
    full_df['amplitude'] = np.abs(full_df['value'])
    full_df['power']     = full_df['amplitude'] ** 2
    full_df['real']      = np.real(full_df['value'])
    full_df['imag']      = np.imag(full_df['value'])
    full_df['harmonic']  = ("l=" + full_df['degree'].astype(str) +
                             " m=" + full_df['order'].astype(str))

    degrees       = np.arange(lmax + 1)
    total_power   = np.zeros(lmax + 1, dtype=np.float64)
    max_amplitude = np.zeros(lmax + 1, dtype=np.float64)

    for l in degrees:
        mask          = full_df['degree'] == l
        amps          = full_df.loc[mask, 'amplitude'].values
        total_power[l]   = np.sum(amps ** 2)
        max_amplitude[l] = np.max(amps) if amps.size > 0 else 0.0

    spectrum_df = pd.DataFrame({
        'degree':         degrees,
        'total_power':    total_power,
        'max_amplitude':  max_amplitude,
        'total_amplitude': np.sqrt(total_power)
    })

    if output_path:
        base, ext = os.path.splitext(output_path)
        if ext.lower() == '.xlsx':
            with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
                full_df.to_excel(writer,
                                 sheet_name='Full_Coefficients', index=False)
                spectrum_df.to_excel(writer,
                                     sheet_name='Power_Spectrum', index=False)
            print(f"Saved: {output_path}")
        elif ext.lower() == '.csv':
            full_df.to_csv(f"{base}_full.csv",     index=False)
            spectrum_df.to_csv(f"{base}_spectrum.csv", index=False)
            print(f"Saved: {base}_full.csv, {base}_spectrum.csv")
        else:
            raise ValueError("Only .xlsx or .csv supported")

    return full_df, spectrum_df


def visualize_power_spectrum(spectrum_df, max_degree=None,
                              log_scale=True, filename=None):
    """Visualize the spherical harmonic power spectrum."""
    df = spectrum_df.copy()
    if max_degree is not None:
        df = df[df['degree'] <= max_degree]

    plt.figure(figsize=(12, 6))
    plt.plot(df['degree'], df['total_power'],
             marker='o', linestyle='-',
             color='#2c7bb6', markersize=6, linewidth=2,
             label='Total Power')
    plt.plot(df['degree'], df['max_amplitude'],
             marker='s', linestyle='--',
             color='#d7191c', markersize=5, linewidth=1.5,
             alpha=0.7, label='Max Amplitude')

    plt.xlabel('Degree (l)', fontsize=12, labelpad=10)
    plt.ylabel('Power / Amplitude' +
               (' (log scale)' if log_scale else ''),
               fontsize=12, labelpad=10)
    plt.title('Spherical Harmonic Power Spectrum', fontsize=14, pad=20)
    plt.xticks(np.arange(0, df['degree'].max() + 1,
                         5 if df['degree'].max() > 20 else 2))
    plt.xlim(-0.5, df['degree'].max() + 0.5)

    if log_scale:
        plt.yscale('log')
        plt.grid(True, which="both", ls="--", alpha=0.3)
    else:
        plt.grid(True, axis='y', ls="--", alpha=0.3)

    plt.legend(fontsize=10, frameon=True, loc='upper right')

    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")
    else:
        plt.show()