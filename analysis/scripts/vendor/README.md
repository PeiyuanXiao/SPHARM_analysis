# Vendored JavaScript libraries

These third-party libraries are embedded (inlined) into the interactive HTML
outputs under `analysis/output/html/` so that every HTML file is fully
self-contained and works offline — no internet connection or CDN is required.

| File | Version | Upstream source | License |
|------|---------|-----------------|---------|
| `three.module.min.js` | r128 | https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.module.min.js | MIT |
| `plotly-2.27.0.min.js` | 2.27.0 | https://cdn.plot.ly/plotly-2.27.0.min.js | MIT |

Consumed by:
- `../SPHARM_modules/spharm_interactive_export.py`      (three.js)
- `../SPHARM_modules/spherical_kde_interactive_export.py` (three.js)
- `../SPHARM_modules/spharm_reconstruction.py`          (three.js)
- `../r_alignment/align_svd.R`                          (plotly)
- `../r_alignment/align_lin2024.R`                      (plotly)
