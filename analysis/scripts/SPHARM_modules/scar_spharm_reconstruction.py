"""
scar_kde_pyvista_viz.py
=======================
Renders static PNG images of spherical harmonic reconstructions
of scar direction KDE data using pyvista.

For each specimen, outputs two side-by-side PNG files:
  - Left : KDE density mapped as deformed sphere
  - Right: SH reconstruction mapped as deformed sphere

Output: /project/analysis/data/derived_data/figures/kde_spharm/
"""

import os
import numpy as np
import pandas as pd
import pyshtools as pysh
import pyvista as pv

# Headless rendering
pv.start_xvfb()

# =============================================================================
# Config
# =============================================================================

DATA_DIR   = "/project/analysis/data/derived_data"
OUTPUT_DIR = f"{DATA_DIR}/figures/kde_spharm"
LMAX       = 20
DH_SIZE    = 64
VIZ_RES    = 60


# =============================================================================
# Load
# =============================================================================

def load_data():
    kde_matrix = np.load(f"{DATA_DIR}/kde_matrix.npy")
    G          = np.load(f"{DATA_DIR}/kde_grid.npy")
    meta       = pd.read_csv(f"{DATA_DIR}/kde_metadata.csv")
    spharm_df  = pd.read_csv(f"{DATA_DIR}/spharm_direction.csv")
    print(f"Loaded {kde_matrix.shape[0]} specimens")
    return kde_matrix, G, meta, spharm_df


# =============================================================================
# Sphere mesh
# =============================================================================

def make_sphere_mesh(res=VIZ_RES):
    colat = np.linspace(0,      np.pi,       res,          endpoint=True)
    lon   = np.linspace(0, 2 * np.pi, 2 * res + 1,    endpoint=True)
    LON, COLAT = np.meshgrid(lon, colat)
    x = np.sin(COLAT) * np.cos(LON)
    y = np.sin(COLAT) * np.sin(LON)
    z = np.cos(COLAT)
    return x, y, z, COLAT, LON


# =============================================================================
# Interpolate KDE onto sphere mesh (density as colour only)
# =============================================================================

def kde_to_density_mesh(kde_vec, G, res=VIZ_RES, kappa=50.0):
    x, y, z, COLAT, LON = make_sphere_mesh(res)
    tx = np.sin(COLAT) * np.cos(LON)
    ty = np.sin(COLAT) * np.sin(LON)
    tz = np.cos(COLAT)

    dot = np.clip(
        tx[:, :, None] * G[:, 0] +
        ty[:, :, None] * G[:, 1] +
        tz[:, :, None] * G[:, 2],
        -1, 1,
    )
    w       = np.exp(kappa * dot)
    density = (w * kde_vec).sum(axis=2) / w.sum(axis=2)
    return density, x, y, z


# =============================================================================
# Build deformed sphere from density grid (pyvista style)
# r = density value → sphere inflated where density is high
# =============================================================================

def density_to_pyvista_mesh(density, x, y, z):
    """
    Scale sphere coordinates by density:
        r = density (normalised to 0.3–1.0 range for visual clarity)
    Returns a pyvista StructuredGrid.
    """
    # Normalise density to [0.3, 1.0] so the shape is clearly deformed
    d_min, d_max = density.min(), density.max()
    if d_max > d_min:
        r = 0.3 + 0.7 * (density - d_min) / (d_max - d_min)
    else:
        r = np.ones_like(density)

    xd = (r * np.sin(np.arccos(np.clip(z, -1, 1))) * np.cos(np.arctan2(y, x)))
    yd = (r * np.sin(np.arccos(np.clip(z, -1, 1))) * np.sin(np.arctan2(y, x)))
    zd = r * z

    # Transpose to match pyvista StructuredGrid convention
    grid = pv.StructuredGrid(xd.T, yd.T, zd.T)
    grid["density"] = density.flatten(order="F")
    return grid


# =============================================================================
# SH reconstruction → density grid
# =============================================================================

def spharm_reconstruct(s_row, lmax=LMAX, res=VIZ_RES):
    coeff_cols  = [c for c in s_row.index if c.startswith("coeff_")]
    coeffs_flat = s_row[coeff_cols].values.astype(np.float64)
    n_expected  = (lmax + 1) ** 2 * 2
    coeffs      = coeffs_flat[:n_expected].reshape(2, lmax + 1, lmax + 1)

    clm     = pysh.SHCoeffs.from_array(coeffs)
    sh_grid = clm.expand(grid="DH", lmax_calc=lmax)
    grid_2d = np.clip(sh_grid.to_array(), 0, None)

    # Build grid unit vectors from actual DH grid shape
    nlat, nlon = grid_2d.shape
    colat_dh   = np.linspace(0,      np.pi,   nlat, endpoint=False)
    lon_dh     = np.linspace(0, 2*np.pi, nlon, endpoint=False)
    L, C       = np.meshgrid(lon_dh, colat_dh)

    G_dh = np.column_stack([
        np.sin(C).ravel() * np.cos(L).ravel(),
        np.sin(C).ravel() * np.sin(L).ravel(),
        np.cos(C).ravel(),
    ])
    v = grid_2d.ravel()
    v = np.clip(v, 0, None)
    s = v.sum()
    if s > 0:
        v /= s

    return kde_to_density_mesh(v, G_dh, res=res, kappa=50.0)


# =============================================================================
# Render and save one specimen
# =============================================================================

def render_specimen(sid, typology, kde_mesh, rec_density, rec_x, rec_y, rec_z,
                    kde_density, x, y, z, output_dir):
    """Render KDE and SH reconstruction side by side and save as PNG."""

    kde_grid = density_to_pyvista_mesh(kde_density, x, y, z)
    rec_grid = density_to_pyvista_mesh(rec_density, rec_x, rec_y, rec_z)

    camera_pos = [(3.5, 3.5, 2.0), (0, 0, 0), (0, 0, 1)]

    pl = pv.Plotter(shape=(1, 2), off_screen=True,
                    window_size=[1400, 600])

    # Left: KDE
    pl.subplot(0, 0)
    pl.add_mesh(kde_grid,
                scalars="density",
                cmap="hot",
                show_edges=False,
                specular=0.5,
                smooth_shading=True)
    pl.add_title(f"KDE density\n{sid}  |  {typology}", font_size=10)
    pl.add_axes(box_args={"color": "gray"})
    pl.camera_position = camera_pos

    # Right: SH reconstruction
    pl.subplot(0, 1)
    pl.add_mesh(rec_grid,
                scalars="density",
                cmap="hot",
                show_edges=False,
                specular=0.5,
                smooth_shading=True)
    pl.add_title(f"SH reconstruction  (lmax={LMAX})\n{sid}", font_size=10)
    pl.add_axes(box_args={"color": "gray"})
    pl.camera_position = camera_pos

    out_path = os.path.join(output_dir, f"{sid}_kde_spharm.png")
    pl.screenshot(out_path, transparent_background=False)
    pl.close()

    return out_path


# =============================================================================
# Main
# =============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    kde_matrix, G, meta, spharm_df = load_data()
    all_ids = list(meta["ID"].astype(str))

    print(f"\nRendering {len(all_ids)} specimens → {OUTPUT_DIR}\n")

    for i, sid in enumerate(all_ids):
        print(f"  [{i+1}/{len(all_ids)}] {sid}...", end="  ")
        try:
            typology = meta["Typology"].iloc[i]
            kde_vec  = kde_matrix[i]
            s_row    = spharm_df[spharm_df["ID"].astype(str) == sid].iloc[0]

            kde_density, x, y, z           = kde_to_density_mesh(kde_vec, G)
            rec_density, rec_x, rec_y, rec_z = spharm_reconstruct(s_row)

            out = render_specimen(
                sid, typology,
                None,
                rec_density, rec_x, rec_y, rec_z,
                kde_density, x, y, z,
                OUTPUT_DIR,
            )
            print(f"✓  →  {os.path.basename(out)}")

        except Exception as e:
            print(f"✗  {e}")

    print(f"\nDone. All figures saved to:\n  {OUTPUT_DIR}")


if __name__ == "__main__":
    main()