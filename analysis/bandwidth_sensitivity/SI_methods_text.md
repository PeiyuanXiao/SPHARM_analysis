## Supplementary Methods — sensitivity of SP-SPHARM to the KDE bandwidth

The SP-SPHARM descriptors depend on one free smoothing parameter: the bandwidth
*h* of the von Mises–Fisher spherical kernel density estimate (vMF-sKDE) applied
to the scar-orientation vectors before the spherical-harmonic expansion, where
the kernel concentration is κ = 1/*h*². The main analysis fixes *h* = 0.35
(κ ≈ 8.16). To verify that the downstream results do not depend on this choice,
we performed a one-parameter sensitivity analysis over *h* ∈ {0.20, 0.25, 0.30,
0.35, 0.40, 0.45, 0.50} (κ ≈ 25.0 down to 4.0).

For each *h* we recomputed the SP-SPHARM power spectra with the same code used in
the main pipeline, changing only the bandwidth and holding everything else fixed:
the SVD-based alignment of scar vectors, the 72 × 36 spherical evaluation grid,
the 64 × 128 Driscoll–Healy grid, the maximum degree l_max = 20, the
amplitude-normalisation of the per-degree power spectrum, and all random seeds.
Because morphology (M-SPHARM) is derived from the surface meshes and is
independent of *h*, it was reused unchanged. As an internal control, the spectra
recomputed at *h* = 0.35 reproduce the cached production values, anchoring the
sweep to the published pipeline.

At each *h* we re-ran the full battery of downstream analyses using the same
estimators as the main study: (i) the per-degree order-selection diagnostics
(cumulative power and cross-specimen coefficient of variation) for the ideal,
experimental, and archaeological assemblages; (ii) the ideal-core standardised
Euclidean-distance structure, summarised by its correlation with the *h* = 0.35
distance matrix and by the persistence of diagnostic separations; (iii) the
experimental-core PERMANOVA by core type (pseudo-F, R², and the Holm-adjusted
pairwise "resolution profile"); and (iv) the morphology–technology decoupling
tests — the global Mantel correlation and the RV coefficient (with permutation
tests) for the experimental and archaeological assemblages, plus the
archaeological scar~core-type PERMANOVA. We then assessed whether each
conclusion — early spectral truncation at l = 1–6, ideal-core discriminability,
experimental core-type separation, and morphology–technology decoupling —
remained qualitatively unchanged across the bandwidth range. Scripts, the per-*h*
metric table, and the diagnostic figures are provided in the research compendium
(`analysis/bandwidth_sensitivity/`).
