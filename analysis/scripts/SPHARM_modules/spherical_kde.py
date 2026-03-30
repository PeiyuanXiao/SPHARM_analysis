import numpy as np
import pandas as pd
 
 
# =============================================================================
# Step 1: Sphere grid
# =============================================================================
 
def make_sphere_grid(n_bearing: int = 72, n_plunge: int = 36) -> np.ndarray:
    """
    Build a uniform evaluation grid on the unit sphere.

    Parameters
    ----------
    n_bearing : int
        Number of azimuthal divisions (longitude), default 72 → 5° steps.
    n_plunge : int
        Number of elevation divisions (latitude), default 36 → 5° steps.
        The grid covers the full sphere from −90° to +90° (poles included).
        A small epsilon is applied only at the exact poles to avoid
        coordinate singularities in downstream spherical-coordinate
        conversions, while keeping effective angular coverage at > 99.9%.

    Returns
    -------
    G : np.ndarray, shape (n_bearing * n_plunge, 3)
        Unit vectors for each grid point (x, y, z).
    """
    # Include poles with a tiny epsilon guard (0.01°) to avoid
    # sin(colat)=0 singularities in spherical-coordinate conversions.
    eps     = np.deg2rad(0.01)
    bearing = np.linspace(0, 2 * np.pi, n_bearing, endpoint=False)
    plunge  = np.linspace(-np.pi / 2 + eps, np.pi / 2 - eps, n_plunge)

    b, p = np.meshgrid(bearing, plunge)
    x = np.cos(p) * np.cos(b)
    y = np.cos(p) * np.sin(b)
    z = np.sin(p)

    G = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    return G
 
 
# =============================================================================
# Step 2: von Mises-Fisher KDE
# =============================================================================
 
def fit_vmf_kde(
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    G: np.ndarray,
    bandwidth: float = 0.35,
) -> np.ndarray:
    """
    Fit a von Mises-Fisher kernel density estimate for one specimen.
 
    For each grid point g in G, the density is:
        density(g) = mean_i [ exp(kappa * dot(g, x_i)) ]
    where kappa = 1 / bandwidth^2.
 
    Parameters
    ----------
    ux, uy, uz : array-like, shape (n_scars,)
        Unit direction vectors of flaking scars (already normalised).
    G : np.ndarray, shape (n_grid, 3)
        Evaluation grid from make_sphere_grid().
    bandwidth : float
        Smoothing bandwidth. Smaller → sharper peaks.
        Typical range: 0.2 (sharp) – 0.5 (smooth).
 
    Returns
    -------
    density : np.ndarray, shape (n_grid,)
        Normalised density values summing to 1.
    """
    X     = np.column_stack([ux, uy, uz]).astype(np.float64)  # (n_scars, 3)
    kappa = 1.0 / bandwidth ** 2
 
    dot_mat = G @ X.T                                          # (n_grid, n_scars)
    density = np.mean(np.exp(kappa * dot_mat), axis=1)
 
    total = density.sum()
    if total == 0:
        raise ValueError("KDE density sums to zero — check input vectors.")
    return density / total
 
 
# =============================================================================
# Step 3: Batch processing
# =============================================================================
 
def batch_spherical_kde(
    directions_df: pd.DataFrame,
    bandwidth: float = 0.35,
    n_bearing: int = 72,
    n_plunge: int = 36,
    id_col: str = "ID",
    ux_col: str = "ux",
    uy_col: str = "uy",
    uz_col: str = "uz",
    typology_col: str = "Typology",
    verbose: bool = True,
) -> dict:
    """
    Run spherical KDE for every specimen in a DataFrame.
 
    Parameters
    ----------
    directions_df : pd.DataFrame
        Must contain columns for ID, ux, uy, uz, and optionally Typology.
    bandwidth : float
        vMF bandwidth (passed to fit_vmf_kde).
    n_bearing, n_plunge : int
        Grid resolution (passed to make_sphere_grid).
    id_col, ux_col, uy_col, uz_col, typology_col : str
        Column name overrides.
    verbose : bool
        Print progress if True.
 
    Returns
    -------
    result : dict with keys
        'kde_matrix'   : np.ndarray (n_specimens, n_grid)
        'G'            : np.ndarray (n_grid, 3)  sphere grid unit vectors
        'ids'          : list of specimen IDs
        'typologies'   : list of typology labels
        'n_scars'      : list of scar counts per specimen
    """
    G       = make_sphere_grid(n_bearing, n_plunge)
    n_grid  = len(G)
    all_ids = directions_df[id_col].unique()
    n_cores = len(all_ids)
 
    kde_matrix   = np.full((n_cores, n_grid), np.nan)
    typologies   = []
    n_scars      = []
 
    if verbose:
        print(f"\nFitting spherical KDE for {n_cores} specimens "
              f"(bandwidth={bandwidth}, grid={n_bearing}×{n_plunge})...")
 
    for i, id_i in enumerate(all_ids):
        df_i = directions_df[directions_df[id_col] == id_i]
 
        kde_matrix[i] = fit_vmf_kde(
            df_i[ux_col].values,
            df_i[uy_col].values,
            df_i[uz_col].values,
            G,
            bandwidth,
        )
 
        typ = df_i[typology_col].iloc[0] if typology_col in df_i.columns else "unknown"
        typologies.append(typ)
        n_scars.append(len(df_i))
 
        if verbose:
            print(f"  [{i+1:>2}/{n_cores}]  ID={str(id_i):<30}  "
                  f"Type={str(typ):<25}  n_scars={len(df_i)}")
 
    if verbose:
        print("Done.\n")
 
    return {
        "kde_matrix": kde_matrix,
        "G":          G,
        "ids":        list(all_ids),
        "typologies": typologies,
        "n_scars":    n_scars,
    }
