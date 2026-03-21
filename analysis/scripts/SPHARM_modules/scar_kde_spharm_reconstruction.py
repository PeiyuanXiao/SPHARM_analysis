"""
scar_kde_spharm_reconstruction.py
==================================
Exports interactive HTML: KDE sphere + SH reconstruction sphere.
Raw numpy arrays are written directly as JS — avoids plotly serialisation issues.
"""

import os
import json
import numpy as np
import pandas as pd
import pyshtools as pysh

DATA_DIR    = "/project/analysis/data/derived_data"
OUTPUT_HTML = f"{DATA_DIR}/kde_spharm_visualisation.html"
LMAX        = 20
DH_SIZE     = 64
VIZ_RES     = 50


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
    colat = np.linspace(0,      np.pi,       res,         endpoint=True)
    lon   = np.linspace(0, 2 * np.pi, 2 * res + 1,  endpoint=True)
    LON, COLAT = np.meshgrid(lon, colat)
    x = np.sin(COLAT) * np.cos(LON)
    y = np.sin(COLAT) * np.sin(LON)
    z = np.cos(COLAT)
    return x, y, z, COLAT, LON


# =============================================================================
# Interpolate density onto mesh
# =============================================================================

def interp_to_mesh(kde_vec, G, res=VIZ_RES, kappa=50.0):
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
# SH reconstruction
# =============================================================================

def spharm_reconstruct(s_row, lmax=LMAX, res=VIZ_RES):
    coeff_cols  = [c for c in s_row.index if c.startswith("coeff_")]
    coeffs_flat = s_row[coeff_cols].values.astype(np.float64)
    n_expected  = (lmax + 1) ** 2 * 2
    coeffs      = coeffs_flat[:n_expected].reshape(2, lmax + 1, lmax + 1)

    clm     = pysh.SHCoeffs.from_array(coeffs)
    sh_grid = clm.expand(grid="DH", lmax_calc=lmax)
    grid_2d = np.clip(sh_grid.to_array(), 0, None)

    nlat, nlon = grid_2d.shape
    colat_dh   = np.linspace(0,      np.pi,   nlat, endpoint=False)
    lon_dh     = np.linspace(0, 2*np.pi, nlon, endpoint=False)
    L, C       = np.meshgrid(lon_dh, colat_dh)

    G_dh = np.column_stack([
        np.sin(C).ravel() * np.cos(L).ravel(),
        np.sin(C).ravel() * np.sin(L).ravel(),
        np.cos(C).ravel(),
    ])
    v    = grid_2d.ravel()
    v    = np.clip(v, 0, None)
    s    = v.sum()
    if s > 0:
        v /= s

    return interp_to_mesh(v, G_dh, res=res, kappa=50.0)


# =============================================================================
# Serialise 2D array → compact JS list-of-lists
# =============================================================================

def arr_to_js(arr: np.ndarray) -> str:
    return json.dumps(arr.tolist())


# =============================================================================
# Export HTML
# =============================================================================

def export_html():
    kde_matrix, G, meta, spharm_df = load_data()
    all_ids = list(meta["ID"].astype(str))

    specimens = {}   # sid -> {x,y,z,kde,rec}

    for i, sid in enumerate(all_ids):
        print(f"  [{i+1}/{len(all_ids)}] {sid}...", end="  ")
        try:
            kde_vec  = kde_matrix[i]
            s_row    = spharm_df[spharm_df["ID"].astype(str) == sid].iloc[0]
            typology = meta["Typology"].iloc[i]

            kde_d, x, y, z        = interp_to_mesh(kde_vec, G)
            rec_d, xr, yr, zr     = spharm_reconstruct(s_row)

            specimens[sid] = {
                "typology" : typology,
                "x"  : arr_to_js(x),
                "y"  : arr_to_js(y),
                "z"  : arr_to_js(z),
                "kde": arr_to_js(kde_d),
                "xr" : arr_to_js(xr),
                "yr" : arr_to_js(yr),
                "zr" : arr_to_js(zr),
                "rec": arr_to_js(rec_d),
            }
            print("✓")
        except Exception as e:
            print(f"✗ {e}")

    # Build JS data block
    js_lines = ["var SPECIMENS = {};"]
    for sid, d in specimens.items():
        js_lines.append(
            f'SPECIMENS[{json.dumps(sid)}] = {{'
            f'"typology":{json.dumps(d["typology"])},'
            f'"x":{d["x"]},"y":{d["y"]},"z":{d["z"]},'
            f'"kde":{d["kde"]},'
            f'"xr":{d["xr"]},"yr":{d["yr"]},"zr":{d["zr"]},'
            f'"rec":{d["rec"]}'
            f'}};'
        )
    js_data = "\n".join(js_lines)

    options_html = "".join(
        f'<option value="{sid}">{sid}</option>'
        for sid in specimens
    )
    ids_json = json.dumps(list(specimens.keys()))

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Spherical KDE &amp; SH Reconstruction</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
  <style>
    body   {{ font-family: sans-serif; background:#f5f7fa; margin:20px; }}
    h2     {{ text-align:center; color:#333; margin-bottom:8px; }}
    .ctrl  {{ text-align:center; margin-bottom:16px; }}
    label  {{ font-size:14px; color:#555; }}
    select {{ font-size:14px; padding:4px 10px; margin-left:8px; }}
    .row   {{ display:flex; justify-content:center; gap:20px; flex-wrap:wrap; }}
    .box   {{ background:white; border-radius:8px;
               box-shadow:0 1px 4px rgba(0,0,0,.1); padding:8px; }}
  </style>
</head>
<body>
  <h2>Spherical KDE &amp; Spherical Harmonic Reconstruction</h2>
  <div class="ctrl">
    <label>Select specimen:</label>
    <select id="sel">{options_html}</select>
  </div>
  <div class="row">
    <div class="box"><div id="plotKDE" style="width:520px;height:500px;"></div></div>
    <div class="box"><div id="plotRec" style="width:520px;height:500px;"></div></div>
  </div>

  <script>
    {js_data}

    var allIds = {ids_json};

    var CAM    = {{eye:{{x:1.5, y:1.5, z:0.9}}}};
    var SCENE  = {{
      xaxis:{{visible:false}},
      yaxis:{{visible:false}},
      zaxis:{{visible:false}},
      aspectmode:"cube",
      bgcolor:"white",
      camera: CAM
    }};

    function render(id) {{
      var d = SPECIMENS[id];
      if (!d) return;

      var traceKDE = {{
        type:"surface",
        x:d.x, y:d.y, z:d.z,
        surfacecolor:d.kde,
        colorscale:"Hot",
        showscale:true,
        colorbar:{{thickness:14, len:0.6}}
      }};
      var traceRec = {{
        type:"surface",
        x:d.xr, y:d.yr, z:d.zr,
        surfacecolor:d.rec,
        colorscale:"Hot",
        showscale:true,
        colorbar:{{thickness:14, len:0.6}}
      }};

      var layoutKDE = {{
        title:{{text:"<b>KDE density</b><br>" + id + " | " + d.typology,
                x:0.5, font:{{size:13}}}},
        scene: SCENE,
        margin:{{t:70,b:10,l:10,r:10}},
        paper_bgcolor:"white"
      }};
      var layoutRec = {{
        title:{{text:"<b>SH reconstruction</b> (lmax={LMAX})<br>" + id,
                x:0.5, font:{{size:13}}}},
        scene: SCENE,
        margin:{{t:70,b:10,l:10,r:10}},
        paper_bgcolor:"white"
      }};

      Plotly.react("plotKDE", [traceKDE], layoutKDE);
      Plotly.react("plotRec", [traceRec], layoutRec);
    }}

    render(allIds[0]);
    document.getElementById("sel").addEventListener("change", function() {{
      render(this.value);
    }});
  </script>
</body>
</html>"""

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nSaved: {OUTPUT_HTML}")
    print(f"Total: {len(specimens)} specimens")


if __name__ == "__main__":
    export_html()