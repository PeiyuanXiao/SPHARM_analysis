import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import re
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform


# ============================================================
# Variance analysis
# ============================================================

def analyze_variance(all_stats, filenames=None, output_path=None):
    """
    Compute and visualize per-degree variance of power spectra.

    Parameters
    ----------
    all_stats : np.ndarray, shape (N, n_degrees)
    filenames : list of str, optional
    output_path : str, optional
        If provided, saves variance CSV to this path.

    Returns
    -------
    var_df : pd.DataFrame
        Columns: degree, variance
    """
    all_stats = np.array(all_stats)
    variances = np.var(all_stats, axis=0)

    print("\n==== Variance Analysis ====")
    print(f"Samples:       {all_stats.shape[0]}")
    print(f"Max variance:  {np.max(variances):.4f} "
          f"(degree {np.argmax(variances)})")
    print(f"Min variance:  {np.min(variances):.4f} "
          f"(degree {np.argmin(variances)})")
    print(f"Mean variance: {np.mean(variances):.4f}")

    sorted_idx = np.argsort(variances)[::-1]
    print("\nTop 5 degrees by variance:")
    for i, idx in enumerate(sorted_idx[:5]):
        print(f"  Rank {i+1}: degree {idx} → {variances[idx]:.4f}")

    plt.figure(figsize=(10, 4))
    plt.plot(variances, marker='o', linewidth=1.5, markersize=4)
    plt.title('Variance per Degree')
    plt.xlabel('Spherical Harmonic Degree')
    plt.ylabel('Variance')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()

    var_df = pd.DataFrame({
        'degree':   np.arange(len(variances)),
        'variance': variances
    })

    if output_path:
        var_df.to_csv(output_path, index=False)
        print(f"Saved variance analysis: {output_path}")

    return var_df


# ============================================================
# PCA
# ============================================================

def analyze_pca(all_stats, filenames=None, n_components=4,
                output_path=None):
    """
    PCA on power spectra with loading matrix and scatter plots.

    Parameters
    ----------
    all_stats : np.ndarray, shape (N, n_degrees)
    filenames : list of str, optional
    n_components : int
    output_path : str, optional
        If provided, saves PCA scores CSV.

    Returns
    -------
    pca_result : np.ndarray, shape (N, n_components)
    pca : sklearn PCA object
    """
    from sklearn.decomposition import PCA

    pca        = PCA(n_components=n_components)
    pca_result = pca.fit_transform(all_stats)

    if output_path and filenames is not None:
        cols = {f'PC{i+1}': pca_result[:, i]
                for i in range(n_components)}
        pd.DataFrame({'filename': filenames, **cols}).to_csv(
            output_path, index=False)
        print(f"Saved PCA scores: {output_path}")

    print("\n==== PCA Loading Matrix ====")
    for i in range(n_components):
        sorted_idx = np.argsort(np.abs(pca.components_[i]))[::-1]
        print(f"\nPC{i+1} ({pca.explained_variance_ratio_[i]:.1%}):")
        for deg in sorted_idx[:5]:
            print(f"  degree {deg}: {pca.components_[i, deg]:.3f}")

    plot_pca_components(pca_result, pca, filenames,
                        x_component=0, y_component=1)
    plot_pca_components(pca_result, pca, filenames,
                        x_component=2, y_component=3)

    return pca_result, pca


def plot_pca_components(pca_result, pca, filenames=None,
                        x_component=0, y_component=1):
    """Scatter plot of two PCA components."""
    plt.figure(figsize=(10, 8))
    plt.scatter(pca_result[:, x_component],
                pca_result[:, y_component],
                alpha=0.7, edgecolors='w', s=40)

    if filenames is not None:
        for i, (x, y) in enumerate(
            zip(pca_result[:, x_component],
                pca_result[:, y_component])
        ):
            plt.text(x, y,
                     os.path.splitext(filenames[i])[0],
                     fontsize=7, alpha=0.8,
                     rotation=20, ha='left', va='center')

    explained = pca.explained_variance_ratio_
    plt.xlabel(f'PC{x_component+1} ({explained[x_component]:.1%})')
    plt.ylabel(f'PC{y_component+1} ({explained[y_component]:.1%})')
    plt.title(f'PCA Components {x_component+1}–{y_component+1}')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


# ============================================================
# Distance matrix
# ============================================================

def calculate_power_distance(all_stats, filenames, output_dir, lmax=20):
    """
    Compute pairwise power spectrum distances and generate a heatmap.

    Uses row normalization (each spectrum sums to 1) to focus on
    relative shape of the power distribution rather than absolute scale.

    Parameters
    ----------
    all_stats : np.ndarray, shape (N, n_degrees)
    filenames : list of str
    output_dir : str
    lmax : int
    """
    power_spectra = all_stats[:, 1:lmax+1]

    # Row normalization
    row_sum    = power_spectra.sum(axis=1, keepdims=True)
    power_norm = power_spectra / np.where(row_sum > 0, row_sum, 1)

    distance_matrix = squareform(pdist(power_norm, metric='cosine'))

    plt.figure(figsize=(12, 10))
    sns.heatmap(distance_matrix,
                annot=len(filenames) <= 30,
                fmt=".2f",
                cmap="viridis",
                xticklabels=filenames,
                yticklabels=filenames)
    plt.title(f"Power Spectrum Distance Heatmap (l=1~{lmax})")
    plt.tight_layout()

    path = os.path.join(output_dir,
                        f"power_distance_heatmap_lmax{lmax}.png")
    plt.savefig(path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Saved heatmap: {path}")

    return distance_matrix


# ============================================================
# UMAP
# ============================================================

def analyze_umap(all_stats, filenames, output_dir, lmax=20,
                 save_plot=True, save_csv=True):
    """
    UMAP dimensionality reduction on power spectra.

    Parameters
    ----------
    all_stats : np.ndarray, shape (N, n_degrees)
    filenames : list of str
    output_dir : str
    lmax : int
    save_plot : bool
    save_csv : bool

    Returns
    -------
    umap_df : pd.DataFrame
        Columns: x, y, filename, category
    """
    power_spectra = all_stats[:, 1:lmax+1]
    n_samples     = power_spectra.shape[0]
    n_neighbors   = min(6, n_samples - 1)

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.06,
        metric='cosine',
        random_state=42
    )
    umap_result = reducer.fit_transform(power_spectra)

    def get_category(label):
        for cat in ["Multifacial", "Subspheroid", "Spheroid", "Polyhedron"]:
            if cat in label:
                return cat
        return "idealmodel"

    clean_names = [
        re.sub(r'-[^-]+$', '', re.sub(r'^\d+_', '', f))
        for f in filenames
    ]

    umap_df = pd.DataFrame({
        'x':        umap_result[:, 0],
        'y':        umap_result[:, 1],
        'filename': clean_names,
        'category': [get_category(f) for f in filenames]
    })

    if save_csv:
        csv_path = os.path.join(output_dir, f"umap_lmax{lmax}.csv")
        umap_df.to_csv(csv_path, index=False)
        print(f"Saved UMAP CSV: {csv_path}")

    if save_plot:
        plt.figure(figsize=(15, 12))
        for i, (x, y) in enumerate(umap_result):
            plt.scatter(x, y, color='black', s=100, alpha=0.8)
            plt.text(x + 0.02, y + 0.02, clean_names[i],
                     fontsize=9, ha='left', va='bottom')

        plt.title(f"UMAP Projection (l=1~{lmax})", fontsize=14)
        plt.xlabel("UMAP Component 1", fontsize=12)
        plt.ylabel("UMAP Component 2", fontsize=12)
        plt.tight_layout()

        plot_path = os.path.join(output_dir, f"umap_lmax{lmax}.png")
        plt.savefig(plot_path, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved UMAP plot: {plot_path}")

    return umap_df


analyze_umap2 = analyze_umap


# ============================================================
# Entry point for post-batch analysis
# ============================================================

def run_batch_analysis(output_csv, output_dir, lmax=20):
    """
    Post-processing analysis after batch SPHARM computation.
    Called after batch_process completes in SPHARM_main.py.

    Runs: variance analysis + UMAP.

    Parameters
    ----------
    output_csv : str
        Path to SPHARM_results.csv.
    output_dir : str
    lmax : int
    """
    df         = pd.read_csv(output_csv)
    power_cols = [c for c in df.columns if c.startswith("power_degree_")]
    X          = df[power_cols].values

    # Variance analysis
    var_path = os.path.join(output_dir, "variance_per_degree.csv")
    analyze_variance(X,
                     filenames=df["specimen_id"].tolist(),
                     output_path=var_path)

    # UMAP
    if len(df) > 1:
        analyze_umap(
            X,
            df["specimen_id"].astype(str).tolist(),
            output_dir,
            lmax=lmax,
            save_plot=True,
            save_csv=True
        )


# ============================================================
# Entry
# ============================================================

if __name__ == "__main__":
    base_dir   = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.abspath(
        os.path.join(base_dir, "..", "..", "..",
                     "data", "drived_data"))
    output_csv = os.path.join(output_dir, "SPHARM_results.csv")

    os.makedirs(output_dir, exist_ok=True)
    run_batch_analysis(output_csv, output_dir, lmax=20)