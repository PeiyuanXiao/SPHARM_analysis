# ==============================================================================
# spherical_kde_interactive_export.py
#
# Standalone: generate interactive HTML (EXP / IM / SDG) that visualises the
# von Mises-Fisher KDE of flaking-scar direction vectors on a single sphere.
#
#   - Single-sphere viewport, density shown via colour map
#   - Specimen dropdown grouped by Typology
#   - View presets / Info panel
#   - Bandwidth slider (live KDE recompute, kappa = 1/bw^2)
#   - Scar-direction scatter overlay toggle
#   - Graticule toggle
#   - Material selector (Flat / Matte / Glossy / Metal)
#   - Type mean (pool all direction vectors of one Typology, then re-fit KDE)
#   - Colour bar + PNG screenshot export
#   - Scheme B: only direction vectors are sent; vMF KDE is computed in JS
#
# Matches the style of spharm_interactive_export.py (same dark theme / fonts).
# ==============================================================================

import os
import json
import base64
from pathlib import Path
import numpy as np
import pandas as pd

# ── Vendored three.js (r128), inlined as a base64 data: URL so the exported
#    HTML is fully self-contained and works offline (no CDN required) ──────────
_THREE_SRC = None

def _three_src() -> str:
    """three.module.min.js (r128) as a data: URL, read once from the vendored copy."""
    global _THREE_SRC
    if _THREE_SRC is None:
        _p = Path(__file__).resolve().parent.parent / "vendor" / "three.module.min.js"
        _THREE_SRC = "data:text/javascript;base64," + base64.b64encode(_p.read_bytes()).decode("ascii")
    return _THREE_SRC

# ── Config ────────────────────────────────────────────────────────────────────
DIR_CSV = "analysis/data/derived_data/directions_aligned_svd.csv"
OUT_DIR = "analysis/output/html/kde_sphere_interactive"

BANDWIDTH = 0.35  # default vMF bandwidth (matches kde_to_spharm_main.py)

os.makedirs(OUT_DIR, exist_ok=True)


# ── Data loading ────────────────────────────────────────────────────────────────

def load_directions(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df['ID'] = df['ID'].astype(str).str.strip()
    missing = {'ID', 'ux', 'uy', 'uz'} - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")
    if 'Typology' not in df.columns:
        df['Typology'] = 'unknown'
    return df


# ── Specimen records ──────────────────────────────────────────────────────────

def build_specimen_data(df, ids):
    records = []
    total = len(ids)
    for n, sid in enumerate(ids, 1):
        sub = df[df['ID'] == sid]
        if sub.empty:
            print(f"  [{n}/{total}] {sid} - skipped (no data)")
            continue
        # Re-normalise direction vectors as a safeguard
        vecs = sub[['ux', 'uy', 'uz']].values.astype(float)
        norms = np.linalg.norm(vecs, axis=1)
        keep = norms > 1e-9
        vecs = vecs[keep] / norms[keep, None]
        if len(vecs) == 0:
            print(f"  [{n}/{total}] {sid} - skipped (all-zero vectors)")
            continue
        typ = str(sub['Typology'].iloc[0])
        entry = {
            'id':   sid,
            'dirs': np.round(vecs, 6).tolist(),
            'meta': {
                'Typology': typ,
                'n_scars':  str(len(vecs)),
            },
        }
        records.append(entry)
        print(f"  [{n}/{total}] {sid} OK  (n_scars={len(vecs)})")
    return records


# ── HTML generation ─────────────────────────────────────────────────────────────

def generate_html(records, group_name, bandwidth):
    data_json = json.dumps(records, separators=(',', ':'))

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>vMF KDE Sphere — {group_name}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Inter:wght@300;400;500;600&display=swap');
  *,*::before,*::after{{box-sizing:border-box;margin:0;padding:0}}
  :root{{
    --bg:#0a0a0f;--sf:#14141e;--sf2:#1e1e2a;--sf3:#282838;
    --bd:#2a2a3a;--bdh:#4a4a6a;--tx:#e4e4f0;--txd:#7878a0;
    --acM:#6cb4d8;--acS:#d8a06c;--acA:#8cd0f0;--r:8px;
  }}
  body{{font-family:'Inter',system-ui,sans-serif;background:var(--bg);
       color:var(--tx);min-height:100vh;overflow:hidden}}

  .hdr{{padding:10px 20px 8px;border-bottom:1px solid var(--bd);
       background:var(--sf);display:flex;flex-direction:column;gap:6px}}
  .hdr-row{{display:flex;align-items:center;gap:8px;flex-wrap:nowrap;overflow-x:auto}}
  .hdr h1{{font-family:'JetBrains Mono',monospace;font-size:14px;
           font-weight:600;letter-spacing:.6px;color:var(--acM);flex-shrink:0}}
  .hdr h1 span{{color:var(--txd);font-weight:400}}

  select{{
    font-family:'JetBrains Mono',monospace;font-size:11px;
    background:var(--sf2);color:var(--tx);border:1px solid var(--bd);
    border-radius:var(--r);padding:5px 24px 5px 8px;cursor:pointer;outline:none;
    appearance:none;
    background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='10' height='6'%3E%3Cpath d='M0 0l5 6 5-6z' fill='%237878a0'/%3E%3C/svg%3E");
    background-repeat:no-repeat;background-position:right 8px center;
    transition:border-color .15s;
  }}
  select:hover{{border-color:var(--bdh)}}
  select:focus{{border-color:var(--acM)}}
  optgroup{{font-style:normal;color:var(--txd)}}
  #specSel{{max-width:280px}}

  .btn{{
    font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:.3px;
    background:var(--sf2);color:var(--txd);border:1px solid var(--bd);
    border-radius:var(--r);padding:5px 10px;cursor:pointer;
    transition:all .15s;white-space:nowrap;user-select:none;
  }}
  .btn:hover{{border-color:var(--bdh);color:var(--tx)}}
  .btn.active{{background:var(--acM);color:var(--bg);border-color:var(--acM);font-weight:600}}

  .sep{{width:1px;height:20px;background:var(--bd);flex-shrink:0}}

  .deg-wrap{{display:flex;align-items:center;gap:5px}}
  .deg-wrap label{{font-family:'JetBrains Mono',monospace;font-size:11px;color:var(--txd);white-space:nowrap}}
  .deg-wrap input[type=range]{{
    width:120px;height:4px;-webkit-appearance:none;appearance:none;
    background:var(--bd);border-radius:2px;outline:none;cursor:pointer;
  }}
  .deg-wrap input[type=range]::-webkit-slider-thumb{{
    -webkit-appearance:none;width:14px;height:14px;border-radius:50%;
    background:var(--acM);border:2px solid var(--bg);cursor:pointer;
  }}
  .deg-wrap input[type=range]::-moz-range-thumb{{
    width:14px;height:14px;border-radius:50%;
    background:var(--acM);border:2px solid var(--bg);cursor:pointer;
  }}
  #bwVal{{
    font-family:'JetBrains Mono',monospace;font-size:12px;font-weight:600;
    color:var(--acA);min-width:34px;text-align:center;
  }}

  .vp-wrap{{display:flex;width:100%;height:calc(100vh - 82px)}}
  .vp{{flex:1;position:relative;overflow:hidden}}
  .vp-lb{{
    position:absolute;top:14px;left:18px;
    font-family:'JetBrains Mono',monospace;font-size:12px;font-weight:600;
    letter-spacing:1.5px;text-transform:uppercase;
    z-index:10;pointer-events:none;user-select:none;opacity:.9;color:var(--acS);
  }}
  .vp-st{{
    position:absolute;top:14px;right:18px;
    font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--txd);
    z-index:10;pointer-events:none;user-select:none;text-align:right;line-height:1.6;
  }}
  .nd{{position:absolute;inset:0;display:flex;align-items:center;
      justify-content:center;font-size:13px;color:var(--txd);font-style:italic}}
  canvas{{display:block}}

  .info-p{{
    position:absolute;bottom:16px;left:18px;
    background:rgba(20,20,30,.94);border:1px solid var(--bd);
    border-radius:var(--r);padding:12px 16px;font-size:12px;line-height:1.7;
    max-width:320px;max-height:200px;overflow-y:auto;z-index:10;display:none;backdrop-filter:blur(8px);
  }}
  .info-p.vis{{display:block}}
  .info-p .mk{{color:var(--txd);font-family:'JetBrains Mono',monospace;font-size:11px;margin-right:6px}}
  .info-p .mv{{color:var(--tx)}}
  .info-p::-webkit-scrollbar{{width:4px}}
  .info-p::-webkit-scrollbar-track{{background:transparent}}
  .info-p::-webkit-scrollbar-thumb{{background:var(--bd);border-radius:2px}}

  .rebuilding{{
    position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);
    font-family:'JetBrains Mono',monospace;font-size:12px;color:var(--acA);
    z-index:20;pointer-events:none;display:none;
  }}
  .rebuilding.vis{{display:block}}

  .colorbar{{
    position:absolute;bottom:16px;right:18px;height:120px;
    z-index:10;pointer-events:none;display:flex;align-items:center;gap:4px;
  }}
  .colorbar canvas{{border-radius:3px;border:1px solid var(--bd)}}
  .cb-ticks{{
    display:flex;flex-direction:column;justify-content:space-between;height:120px;
    font-family:'JetBrains Mono',monospace;font-size:8px;color:var(--txd);
  }}

  .map2d{{
    position:absolute;bottom:16px;left:50%;transform:translateX(-50%);
    background:rgba(20,20,30,.94);border:1px solid var(--bd);border-radius:var(--r);
    padding:8px 10px 10px;z-index:11;display:none;backdrop-filter:blur(8px);
  }}
  .map2d.vis{{display:block}}
  .map2d-hd{{font-family:'JetBrains Mono',monospace;font-size:10px;color:var(--txd);margin-bottom:6px}}
  .map2d canvas{{display:block;border:1px solid var(--bd);border-radius:3px}}
</style>
</head>
<body>

<div class="hdr">
 <div class="hdr-row">
  <h1>vMF KDE <span>● {group_name}</span></h1>
  <select id="specSel"></select>
  <div class="sep"></div>
  <div class="deg-wrap">
    <label>bw</label>
    <input type="range" id="bwSlider" min="0.10" max="0.80" value="{bandwidth}" step="0.01">
    <span id="bwVal">{bandwidth:.2f}</span>
  </div>
  <div class="sep"></div>
  <select id="matSel" title="Surface material">
    <option value="flat">Flat</option>
    <option value="mosaic">Mosaic</option>
    <option value="matte">Matte</option>
    <option value="glossy">Glossy</option>
    <option value="metal">Metal</option>
  </select>
  <select id="cmapSel" title="Colormap scheme">
    <option value="viridis" selected>Viridis</option>
    <option value="magma">Magma</option>
    <option value="plasma">Plasma</option>
    <option value="inferno">Inferno</option>
  </select>
 </div>
 <div class="hdr-row">
  <button class="btn active" data-view="iso">Iso</button>
  <button class="btn" data-view="top">Top</button>
  <button class="btn" data-view="front">Front</button>
  <div class="sep"></div>
  <button class="btn active" id="btnPts">Scars</button>
  <button class="btn" id="btnGrid">Graticule</button>
  <button class="btn" id="btnInfo">Info</button>
  <button class="btn" id="btnMean">Type Mean</button>
  <div class="sep"></div>
  <button class="btn" id="btnContours">Contours</button>
  <button class="btn" id="btnMap">2D Map</button>
  <div class="sep"></div>
  <button class="btn" id="btnSnap" title="Save screenshot">📷 PNG</button>
 </div>
</div>

<div class="vp-wrap">
  <div class="vp" id="vp">
    <span class="vp-lb">Scar Direction KDE</span>
    <div class="vp-st" id="st"></div>
    <div class="nd" id="nd" style="display:none">No data</div>
    <div class="info-p" id="info"></div>
    <div class="rebuilding" id="rb">Computing KDE...</div>
    <div class="colorbar"><canvas width="16" height="120" id="cbc"></canvas><div class="cb-ticks" id="cbTicks"></div></div>
    <div class="map2d" id="map2d"><div class="map2d-hd">2D density — bearing (x, 0–360°) × plunge (y, +90°…−90°)</div><canvas id="mapc" width="360" height="180"></canvas></div>
  </div>
</div>

<script type="importmap">
{{
  "imports": {{
    "three":"{_three_src()}"
  }}
}}
</script>

<script type="module">
import * as THREE from 'three';

const NLAT=130;   // sphere grid latitude resolution
const NLON=260;   // sphere grid longitude resolution

/* ================================================================
   Colour map: viridis (5-anchor linear interpolation), low -> high density
   ================================================================ */
const PALETTE=[
  [0.267,0.005,0.329],
  [0.231,0.322,0.545],
  [0.128,0.567,0.551],
  [0.369,0.789,0.383],
  [0.993,0.906,0.144],
];
// extra sequential maps (6th-degree polynomial approximations of matplotlib
// colormaps, coefficients after Matt Zucker, shadertoy WlfXRN)
const CMAP_POLY={{
  magma:[[-0.00214,-0.00075,-0.00539],[0.25166,0.67752,2.49403],[8.35372,-3.57772,0.31447],[-27.66873,14.26473,-13.64921],[52.17614,-27.94361,12.94417],[-50.76853,29.04658,4.23415],[18.65571,-11.48977,-5.60196]],
  plasma:[[0.05873,0.02334,0.54334],[2.17651,0.23838,0.75396],[-2.68946,-7.45585,3.11080],[6.13035,42.34619,-28.51885],[-11.10744,-82.66631,60.13985],[10.02307,71.41362,-54.07219],[-3.65871,-22.93153,18.19191]],
  inferno:[[0.00022,0.00165,-0.01948],[0.10651,0.56396,3.93271],[11.60249,-3.97285,-15.94239],[-41.70400,17.43640,44.35415],[77.16294,-33.40236,-81.80731],[-71.31943,32.62606,73.20952],[25.13113,-12.24267,-23.07033]]
}};
function cmap(t){{
  t=Math.max(0,Math.min(1,t));
  if(curCmap==='viridis'){{ // original 5-anchor viridis (default, unchanged)
    const s=t*(PALETTE.length-1);
    const i=Math.min(PALETTE.length-2,Math.floor(s));
    const f=s-i;
    const a=PALETTE[i],b=PALETTE[i+1];
    return [a[0]+(b[0]-a[0])*f,a[1]+(b[1]-a[1])*f,a[2]+(b[2]-a[2])*f];
  }}
  const c=CMAP_POLY[curCmap]||CMAP_POLY.magma;
  const out=[0,0,0];
  for(let j=0;j<3;j++){{
    let v=c[6][j];
    for(let p=5;p>=0;p--)v=v*t+c[p][j];
    out[j]=v<0?0:(v>1?1:v);
  }}
  return out;
}}

const MOSAIC_LEVELS=7;                   // discrete colour bands for Mosaic mode
const CONTOUR_LEVELS=[0.25,0.5,0.75];    // iso-density levels (fraction of range)

// quantise t into bands when Mosaic material is active, else pass through
function applyMode(t){{
  if(curMat!=='mosaic')return t;
  const band=Math.min(MOSAIC_LEVELS-1,Math.floor(t*MOSAIC_LEVELS));
  return (band+0.5)/MOSAIC_LEVELS;
}}

// marching-squares edge connections, keyed by 4-corner case (c0=1,c1=2,c2=4,c3=8)
// edges: e0 top(c0-c1) e1 right(c1-c2) e2 bottom(c2-c3) e3 left(c3-c0)
const MS_SEGS={{
  1:[['e3','e0']],2:[['e0','e1']],3:[['e3','e1']],4:[['e1','e2']],
  5:[['e3','e0'],['e1','e2']],6:[['e0','e2']],7:[['e3','e2']],8:[['e2','e3']],
  9:[['e0','e2']],10:[['e0','e1'],['e2','e3']],11:[['e1','e2']],12:[['e1','e3']],
  13:[['e0','e1']],14:[['e0','e3']],
}};

/* ================================================================
   Unit-sphere grid (vertex position == direction vector g)
   ================================================================ */
function buildSphereGrid(nLat,nLon){{
  const nV=nLat*nLon;
  const pos=new Float32Array(nV*3);
  for(let i=0;i<nLat;i++){{
    const theta=(i/(nLat-1))*Math.PI;          // colatitude 0..pi
    const sT=Math.sin(theta),cT=Math.cos(theta);
    for(let j=0;j<nLon;j++){{
      const phi=(j/(nLon-1))*2*Math.PI;
      const idx=i*nLon+j;
      pos[idx*3]  =sT*Math.cos(phi);
      pos[idx*3+1]=sT*Math.sin(phi);
      pos[idx*3+2]=cT;
    }}
  }}
  const indices=[];
  for(let i=0;i<nLat-1;i++)for(let j=0;j<nLon-1;j++){{
    const a=i*nLon+j,b=a+1,c=a+nLon,d=c+1;
    indices.push(a,c,b,b,c,d);
  }}
  return {{pos,indices,nV}};
}}

/* ================================================================
   vMF KDE: density(g)=mean_i exp(kappa*dot(g,x_i)), then normalise
   ================================================================ */
function computeKDE(pos,nV,dirs,bandwidth){{
  const kappa=1.0/(bandwidth*bandwidth);
  const n=dirs.length;
  const dens=new Float64Array(nV);
  // split direction components
  const dx=new Float64Array(n),dy=new Float64Array(n),dz=new Float64Array(n);
  for(let k=0;k<n;k++){{dx[k]=dirs[k][0];dy[k]=dirs[k][1];dz[k]=dirs[k][2];}}
  let total=0;
  for(let v=0;v<nV;v++){{
    const gx=pos[v*3],gy=pos[v*3+1],gz=pos[v*3+2];
    let s=0;
    for(let k=0;k<n;k++)s+=Math.exp(kappa*(gx*dx[k]+gy*dy[k]+gz*dz[k]));
    const d=s/n;
    dens[v]=d;total+=d;
  }}
  if(total>0)for(let v=0;v<nV;v++)dens[v]/=total;
  return dens;
}}

/* ================================================================
   OrbitControls (same as spharm_interactive_export)
   ================================================================ */
class OC{{
  constructor(cam,el){{
    this.cam=cam;this.el=el;this.target=new THREE.Vector3();
    this.sp=new THREE.Spherical();this.spD=new THREE.Spherical();
    this._sc=1;this._pan=new THREE.Vector3();
    this._rs=new THREE.Vector2();this._ps=new THREE.Vector2();
    this._st=0;this._ptrs=[];this._pp={{}};
    this.dampF=0.08;this.rotSpd=0.8;this.panSpd=0.6;
    const off=new THREE.Vector3().copy(cam.position).sub(this.target);
    this.sp.setFromVector3(off);
    el.addEventListener('pointerdown',e=>this._pd(e));
    el.addEventListener('pointermove',e=>this._pm(e));
    el.addEventListener('pointerup',e=>this._pu(e));
    el.addEventListener('pointercancel',e=>this._pu(e));
    el.addEventListener('wheel',e=>{{e.preventDefault();this._sc*=e.deltaY>0?1.06:0.94;}},{{passive:false}});
    el.addEventListener('contextmenu',e=>e.preventDefault());
  }}
  _pd(e){{this.el.setPointerCapture(e.pointerId);this._ptrs.push(e.pointerId);this._pp[e.pointerId]=new THREE.Vector2(e.clientX,e.clientY);if(this._ptrs.length===1){{this._st=e.button===2?3:1;(this._st===1?this._rs:this._ps).set(e.clientX,e.clientY);}}else if(this._ptrs.length===2){{this._st=4;const a=this._pp[this._ptrs[0]],b=this._pp[this._ptrs[1]];this._pinch0=a.distanceTo(b);}}}}
  _pm(e){{if(!this._pp[e.pointerId])return;this._pp[e.pointerId].set(e.clientX,e.clientY);const h=this.el.clientHeight||1;if(this._st===1){{this.spD.theta-=(e.clientX-this._rs.x)/h*Math.PI*this.rotSpd;this.spD.phi-=(e.clientY-this._rs.y)/h*Math.PI*this.rotSpd;this._rs.set(e.clientX,e.clientY);}}else if(this._st===3){{const off=new THREE.Vector3().copy(this.cam.position).sub(this.target);const d=off.length();const up=new THREE.Vector3().copy(this.cam.up).normalize();const rt=new THREE.Vector3().crossVectors(up,off).normalize();this._pan.add(rt.multiplyScalar(-(e.clientX-this._ps.x)/h*this.panSpd*d));this._pan.add(up.multiplyScalar((e.clientY-this._ps.y)/h*this.panSpd*d));this._ps.set(e.clientX,e.clientY);}}else if(this._st===4&&this._ptrs.length===2){{const a=this._pp[this._ptrs[0]],b=this._pp[this._ptrs[1]];const d=a.distanceTo(b);this._sc*=this._pinch0/d;this._pinch0=d;}}}}
  _pu(e){{try{{this.el.releasePointerCapture(e.pointerId);}}catch{{}}this._ptrs=this._ptrs.filter(id=>id!==e.pointerId);delete this._pp[e.pointerId];if(!this._ptrs.length)this._st=0;}}
  update(){{const off=new THREE.Vector3().copy(this.cam.position).sub(this.target);this.sp.setFromVector3(off);this.sp.theta+=this.spD.theta*this.dampF;this.sp.phi+=this.spD.phi*this.dampF;this.sp.phi=Math.max(0.01,Math.min(Math.PI-0.01,this.sp.phi));this.sp.radius*=this._sc;this.sp.radius=Math.max(0.05,Math.min(80,this.sp.radius));this.target.add(this._pan);off.setFromSpherical(this.sp);this.cam.position.copy(this.target).add(off);this.cam.lookAt(this.target);this.spD.theta*=(1-this.dampF);this.spD.phi*=(1-this.dampF);this._sc=1;this._pan.set(0,0,0);}}
  getS(){{return{{t:this.sp.theta,p:this.sp.phi,r:this.sp.radius,tg:this.target.clone()}};}}
  setS(s){{this.sp.theta=s.t;this.sp.phi=s.p;this.sp.radius=s.r;this.target.copy(s.tg);const off=new THREE.Vector3().setFromSpherical(this.sp);this.cam.position.copy(this.target).add(off);this.cam.lookAt(this.target);this.spD.set(0,0,0);this._sc=1;this._pan.set(0,0,0);}}
}}

/* ================================================================
   Data & State
   ================================================================ */
const DATA={data_json};
let curIdx=0,curBW={bandwidth},curMat='flat',curCmap='viridis';
let showPts=true,showGrid=false,showInfo=false,showMean=false,showContours=false,showMap=false;

const SPHERE=buildSphereGrid(NLAT,NLON);

/* ================================================================
   Viewport
   ================================================================ */
const ctr=document.getElementById('vp');
const scene=new THREE.Scene();scene.background=new THREE.Color(0x0a0a0f);
const cam=new THREE.PerspectiveCamera(32,1,0.01,100);cam.position.set(2.5,1.8,2.5);
const rdr=new THREE.WebGLRenderer({{antialias:true,alpha:false,preserveDrawingBuffer:true,powerPreference:'high-performance'}});
rdr.setPixelRatio(Math.min(window.devicePixelRatio,2));
rdr.outputEncoding=THREE.sRGBEncoding;
ctr.appendChild(rdr.domElement);
const ctrl=new OC(cam,rdr.domElement);

// lights (lit materials only; Flat uses MeshBasicMaterial and ignores them)
scene.add(new THREE.AmbientLight(0xffffff,0.65));
const lt1=new THREE.DirectionalLight(0xffffff,0.7);lt1.position.set(4,6,5);scene.add(lt1);
const lt2=new THREE.DirectionalLight(0xc0c8e0,0.3);lt2.position.set(-4,3,-3);scene.add(lt2);

// material factory (vertex colours always carry KDE density)
function makeSphereMat(key){{
  if(key==='flat'||key==='mosaic')   // unlit, so colour bands read as pure density
    return new THREE.MeshBasicMaterial({{vertexColors:true,side:THREE.DoubleSide}});
  const d={{vertexColors:true,side:THREE.DoubleSide}};
  if(key==='matte'){{d.roughness=0.92;d.metalness=0.0;}}
  else if(key==='glossy'){{d.roughness=0.30;d.metalness=0.10;}}
  else if(key==='metal'){{d.roughness=0.40;d.metalness=0.85;}}
  return new THREE.MeshStandardMaterial(d);
}}

// sphere geometry (fixed; only vertex colours are updated)
const sphereGeom=new THREE.BufferGeometry();
sphereGeom.setAttribute('position',new THREE.BufferAttribute(SPHERE.pos,3));
sphereGeom.setIndex(SPHERE.indices);
sphereGeom.computeVertexNormals();
const sphereColors=new Float32Array(SPHERE.nV*3);
sphereGeom.setAttribute('color',new THREE.BufferAttribute(sphereColors,3));
const sphereMesh=new THREE.Mesh(sphereGeom,makeSphereMat(curMat));
scene.add(sphereMesh);

// graticule
const grat=makeGraticule();grat.visible=showGrid;scene.add(grat);

// scar scatter / iso-density contours / cached density grid
let ptsObj=null,contourObj=null,lastDens=null;

function makeGraticule(){{
  const g=new THREE.Group();
  const mat=new THREE.LineBasicMaterial({{color:0x4a4a6a,transparent:true,opacity:0.45}});
  const R=1.003;
  // parallels (every 15 deg)
  for(let lat=-75;lat<=75;lat+=15){{
    const th=lat*Math.PI/180,pts=[];
    for(let a=0;a<=360;a+=3){{const ph=a*Math.PI/180;pts.push(new THREE.Vector3(R*Math.cos(th)*Math.cos(ph),R*Math.cos(th)*Math.sin(ph),R*Math.sin(th)));}}
    g.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),mat));
  }}
  // meridians (every 15 deg)
  for(let lon=0;lon<360;lon+=15){{
    const ph=lon*Math.PI/180,pts=[];
    for(let a=-90;a<=90;a+=4){{const th=a*Math.PI/180;pts.push(new THREE.Vector3(R*Math.cos(th)*Math.cos(ph),R*Math.cos(th)*Math.sin(ph),R*Math.sin(th)));}}
    g.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts),mat));
  }}
  return g;
}}

function buildPoints(dirs){{
  if(ptsObj){{scene.remove(ptsObj);ptsObj.geometry.dispose();ptsObj.material.dispose();ptsObj=null;}}
  if(!dirs||!dirs.length)return;
  const p=new Float32Array(dirs.length*3);
  for(let k=0;k<dirs.length;k++){{p[k*3]=dirs[k][0]*1.012;p[k*3+1]=dirs[k][1]*1.012;p[k*3+2]=dirs[k][2]*1.012;}}
  const geo=new THREE.BufferGeometry();
  geo.setAttribute('position',new THREE.BufferAttribute(p,3));
  const mat=new THREE.PointsMaterial({{color:0x8cd0f0,size:0.045,sizeAttenuation:true}});
  ptsObj=new THREE.Points(geo,mat);ptsObj.visible=showPts;scene.add(ptsObj);
}}

/* direction vectors to display (single specimen or pooled type mean) */
function currentDirs(){{
  const cur=DATA[curIdx];
  if(!showMean)return cur.dirs;
  const typ=cur.meta&&cur.meta.Typology;
  if(!typ)return cur.dirs;
  const merged=[];
  for(const d of DATA)if(d.meta&&d.meta.Typology===typ)for(const v of d.dirs)merged.push(v);
  return merged.length?merged:cur.dirs;
}}

let dMin=0,dMax=0;
function recolor(dirs){{
  const dens=computeKDE(SPHERE.pos,SPHERE.nV,dirs,curBW);
  let mx=-Infinity,mn=Infinity;
  for(let v=0;v<SPHERE.nV;v++){{if(dens[v]>mx)mx=dens[v];if(dens[v]<mn)mn=dens[v];}}
  dMin=mn;dMax=mx;lastDens=dens;
  const rng=(mx-mn)||1;
  const col=sphereGeom.getAttribute('color');
  for(let v=0;v<SPHERE.nV;v++){{
    const c=cmap(applyMode((dens[v]-mn)/rng));
    col.array[v*3]=c[0];col.array[v*3+1]=c[1];col.array[v*3+2]=c[2];
  }}
  col.needsUpdate=true;
  buildContours(dens);
  drawMap2D();
}}

/* ================================================================
   Iso-density contours on the sphere (marching squares on the grid)
   ================================================================ */
function buildContours(dens){{
  if(contourObj){{scene.remove(contourObj);contourObj.geometry.dispose();contourObj.material.dispose();contourObj=null;}}
  if(!showContours||!dens)return;
  const rng=(dMax-dMin)||1,R=1.006,pos=SPHERE.pos,verts=[];
  const P=idx=>[pos[idx*3],pos[idx*3+1],pos[idx*3+2]];
  function cross(iA,iB,L){{                       // crossing point on edge iA-iB, projected to sphere
    let s=(L-dens[iA])/(dens[iB]-dens[iA]);s=Math.max(0,Math.min(1,s));
    const a=P(iA),b=P(iB);
    let x=a[0]+(b[0]-a[0])*s,y=a[1]+(b[1]-a[1])*s,z=a[2]+(b[2]-a[2])*s;
    const n=Math.hypot(x,y,z)||1;return [x/n*R,y/n*R,z/n*R];
  }}
  for(const t of CONTOUR_LEVELS){{
    const L=dMin+t*rng;
    for(let i=0;i<NLAT-1;i++)for(let j=0;j<NLON-1;j++){{
      const i0=i*NLON+j,i1=i*NLON+j+1,i2=(i+1)*NLON+j+1,i3=(i+1)*NLON+j;
      let cs=0;
      if(dens[i0]>=L)cs|=1;if(dens[i1]>=L)cs|=2;if(dens[i2]>=L)cs|=4;if(dens[i3]>=L)cs|=8;
      const segs=MS_SEGS[cs];if(!segs)continue;
      const e={{}};
      if((dens[i0]>=L)!==(dens[i1]>=L))e.e0=cross(i0,i1,L);
      if((dens[i1]>=L)!==(dens[i2]>=L))e.e1=cross(i1,i2,L);
      if((dens[i2]>=L)!==(dens[i3]>=L))e.e2=cross(i2,i3,L);
      if((dens[i3]>=L)!==(dens[i0]>=L))e.e3=cross(i3,i0,L);
      for(const [p,q] of segs){{const A=e[p],B=e[q];if(A&&B)verts.push(A[0],A[1],A[2],B[0],B[1],B[2]);}}
    }}
  }}
  if(!verts.length)return;
  const geo=new THREE.BufferGeometry();
  geo.setAttribute('position',new THREE.Float32BufferAttribute(verts,3));
  contourObj=new THREE.LineSegments(geo,new THREE.LineBasicMaterial({{color:0xffffff,transparent:true,opacity:0.55}}));
  scene.add(contourObj);
}}

/* ================================================================
   2D equirectangular density map (bearing × plunge)
   ================================================================ */
function drawMap2D(){{
  if(!showMap||!lastDens)return;
  const c=document.getElementById('mapc'),ctx=c.getContext('2d');
  const W=c.width,H=c.height,rng=(dMax-dMin)||1;
  const img=ctx.createImageData(W,H);
  for(let y=0;y<H;y++){{
    const i=Math.min(NLAT-1,Math.round(y/(H-1)*(NLAT-1)));   // row 0 = north pole (plunge +90)
    for(let x=0;x<W;x++){{
      const j=Math.min(NLON-1,Math.round(x/(W-1)*(NLON-1)));
      const col=cmap(applyMode((lastDens[i*NLON+j]-dMin)/rng));
      const o=(y*W+x)*4;
      img.data[o]=col[0]*255;img.data[o+1]=col[1]*255;img.data[o+2]=col[2]*255;img.data[o+3]=255;
    }}
  }}
  ctx.putImageData(img,0,0);
  if(showPts){{                                                // overlay scar directions
    ctx.fillStyle='#8cd0f0';
    for(const v of currentDirs()){{
      let bearing=Math.atan2(v[1],v[0]);if(bearing<0)bearing+=2*Math.PI;
      const plunge=Math.asin(Math.max(-1,Math.min(1,v[2])));
      const x=bearing/(2*Math.PI)*W,yy=(1-(plunge+Math.PI/2)/Math.PI)*H;
      ctx.beginPath();ctx.arc(x,yy,1.6,0,2*Math.PI);ctx.fill();
    }}
  }}
}}

/* ================================================================
   Colour bar
   ================================================================ */
function drawColorbar(){{
  const c=document.getElementById('cbc');const ctx=c.getContext('2d');
  const h=c.height,w=c.width;
  for(let y=0;y<h;y++){{
    const col=cmap(applyMode(1-y/h));
    ctx.fillStyle=`rgb(${{Math.round(col[0]*255)}},${{Math.round(col[1]*255)}},${{Math.round(col[2]*255)}})`;
    ctx.fillRect(0,y,w,1);
  }}
  const rng=(dMax-dMin)||0;                       // tick labels: density at t = 1,.75,.5,.25,0
  document.getElementById('cbTicks').innerHTML=
    [1,0.75,0.5,0.25,0].map(t=>`<span>${{(dMin+t*rng).toExponential(1)}}</span>`).join('');
}}

/* ================================================================
   Load specimen
   ================================================================ */
function loadSpec(idx){{
  curIdx=idx;
  const r=DATA[idx];
  const dirs=currentDirs();
  document.getElementById('nd').style.display=(dirs&&dirs.length)?'none':'flex';
  recolor(dirs);
  buildPoints(dirs);
  drawColorbar();
  updateStats(dirs);
  updateInfo(r,dirs);
}}

function updateStats(dirs){{
  const lbl=showMean?' (mean)':'';
  document.getElementById('st').innerHTML=
    `n_scars: ${{dirs.length}}${{lbl}}<br>bw: ${{curBW.toFixed(2)}} (κ=${{(1/(curBW*curBW)).toFixed(1)}})<br>grid: ${{NLAT}}×${{NLON}}`;
}}

function updateInfo(r,dirs){{
  let extra='';
  if(showMean){{
    const typ=r.meta&&r.meta.Typology;
    if(typ){{const nSp=DATA.filter(d=>d.meta&&d.meta.Typology===typ).length;
      extra=`<br><span class="mk">Mode:</span><span class="mv">Type mean (${{typ}}, ${{nSp}} specimens, ${{dirs.length}} scars)</span>`;}}
  }}
  const idLine=`<span class="mk">ID:</span><span class="mv">${{r.id}}</span>`;
  const metaHtml=Object.entries(r.meta||{{}}).map(([k,v])=>`<span class="mk">${{k}}:</span><span class="mv">${{v}}</span>`).join('<br>');
  document.getElementById('info').innerHTML=idLine+(metaHtml?'<br>'+metaHtml:'')+extra;
}}

/* ================================================================
   Dropdown (grouped by Typology)
   ================================================================ */
const sel=document.getElementById('specSel');
(function(){{
  const groups={{}};
  DATA.forEach((d,i)=>{{const typ=(d.meta&&d.meta.Typology)||'Other';if(!groups[typ])groups[typ]=[];groups[typ].push({{i,id:d.id}});}});
  const keys=Object.keys(groups).sort();
  if(keys.length===1){{
    DATA.forEach((d,i)=>{{const o=document.createElement('option');o.value=i;o.textContent=d.id;sel.appendChild(o);}});
  }}else{{
    keys.forEach(k=>{{const og=document.createElement('optgroup');og.label=k;
      groups[k].forEach(x=>{{const o=document.createElement('option');o.value=x.i;o.textContent=x.id;og.appendChild(o);}});
      sel.appendChild(og);}});
  }}
}})();
sel.addEventListener('change',()=>loadSpec(parseInt(sel.value)));

/* ================================================================
   Bandwidth slider (debounced recompute)
   ================================================================ */
const bwSlider=document.getElementById('bwSlider');
const bwVal=document.getElementById('bwVal');
let bwTimer=null;
bwSlider.addEventListener('input',()=>{{
  curBW=parseFloat(bwSlider.value);bwVal.textContent=curBW.toFixed(2);
  document.getElementById('rb').classList.add('vis');
  clearTimeout(bwTimer);
  bwTimer=setTimeout(()=>{{
    const dirs=currentDirs();
    recolor(dirs);drawColorbar();updateStats(dirs);
    document.getElementById('rb').classList.remove('vis');
  }},60);
}});

/* ================================================================
   View presets
   ================================================================ */
const VIEWS={{iso:{{t:Math.PI/4,p:Math.PI/3}},top:{{t:0,p:0.01}},front:{{t:Math.PI/2,p:Math.PI/2}}}};
document.querySelectorAll('[data-view]').forEach(b=>{{
  b.addEventListener('click',()=>{{
    document.querySelectorAll('[data-view]').forEach(x=>x.classList.remove('active'));
    b.classList.add('active');const v=VIEWS[b.dataset.view];
    const s=ctrl.getS();s.t=v.t;s.p=v.p;ctrl.setS(s);
  }});
}});

/* ================================================================
   Toggles
   ================================================================ */
const bPts=document.getElementById('btnPts');
bPts.addEventListener('click',()=>{{showPts=!showPts;bPts.classList.toggle('active',showPts);if(ptsObj)ptsObj.visible=showPts;drawMap2D();}});

const bGrid=document.getElementById('btnGrid');
bGrid.addEventListener('click',()=>{{showGrid=!showGrid;bGrid.classList.toggle('active',showGrid);grat.visible=showGrid;}});

const bContours=document.getElementById('btnContours');
bContours.addEventListener('click',()=>{{showContours=!showContours;bContours.classList.toggle('active',showContours);buildContours(lastDens);}});

const bMap=document.getElementById('btnMap');
bMap.addEventListener('click',()=>{{showMap=!showMap;bMap.classList.toggle('active',showMap);document.getElementById('map2d').classList.toggle('vis',showMap);drawMap2D();}});

const bInfo=document.getElementById('btnInfo');
bInfo.addEventListener('click',()=>{{showInfo=!showInfo;bInfo.classList.toggle('active',showInfo);document.getElementById('info').classList.toggle('vis',showInfo);}});

const bMean=document.getElementById('btnMean');
bMean.addEventListener('click',()=>{{showMean=!showMean;bMean.classList.toggle('active',showMean);loadSpec(curIdx);}});

/* ================================================================
   Material selector (swap material, keep geometry and density colours)
   ================================================================ */
document.getElementById('matSel').addEventListener('change',e=>{{
  curMat=e.target.value;
  const old=sphereMesh.material;
  sphereMesh.material=makeSphereMat(curMat);
  old.dispose();
  recolor(currentDirs());drawColorbar();   // Mosaic changes the colour banding
}});

/* ================================================================
   Colormap scheme selector
   ================================================================ */
document.getElementById('cmapSel').addEventListener('change',e=>{{
  curCmap=e.target.value;
  recolor(currentDirs());drawColorbar();
}});

/* ================================================================
   Export PNG
   ================================================================ */
document.getElementById('btnSnap').addEventListener('click',()=>{{
  rdr.render(scene,cam);
  const data=rdr.domElement.toDataURL('image/png');
  const img=new Image();
  img.onload=()=>{{
    const c=document.createElement('canvas');c.width=img.width;c.height=img.height;
    const ctx=c.getContext('2d');ctx.fillStyle='#0a0a0f';ctx.fillRect(0,0,c.width,c.height);
    ctx.drawImage(img,0,0);
    ctx.font='bold 16px JetBrains Mono, monospace';ctx.fillStyle='#d8a06c';
    ctx.fillText('Scar Direction KDE',20,30);
    ctx.font='12px JetBrains Mono, monospace';ctx.fillStyle='#7878a0';
    const info=DATA[curIdx].id+(showMean?' (type mean)':'')+' | bw='+curBW.toFixed(2);
    ctx.fillText(info,20,c.height-15);
    c.toBlob(blob=>{{
      if(!blob)return;
      const url=URL.createObjectURL(blob);const a=document.createElement('a');
      a.href=url;a.download='kde_'+DATA[curIdx].id.replace(/[\\s\\/]/g,'_')+'_bw'+curBW.toFixed(2)+'.png';
      document.body.appendChild(a);a.click();document.body.removeChild(a);
      setTimeout(()=>URL.revokeObjectURL(url),1000);
    }},'image/png');
  }};
  img.src=data;
}});

/* ================================================================
   Resize & render
   ================================================================ */
function onResize(){{const w=ctr.clientWidth,h=ctr.clientHeight;if(w<1||h<1)return;cam.aspect=w/h;cam.updateProjectionMatrix();rdr.setSize(w,h);}}
window.addEventListener('resize',onResize);onResize();

function animate(){{requestAnimationFrame(animate);ctrl.update();rdr.render(scene,cam);}}

loadSpec(0);
animate();
</script>
</body>
</html>'''
    return html


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("vMF KDE sphere - interactive HTML export")
    print("=" * 60)

    print("\nLoading direction vectors...")
    df = load_directions(DIR_CSV)
    print(f"  {df['ID'].nunique()} specimens, {len(df)} direction vectors")

    all_ids = sorted(df['ID'].unique().tolist())
    groups = {
        'EXP': sorted([i for i in all_ids if i.startswith('EXP')]),
        'IM':  sorted([i for i in all_ids if i.startswith('IM_')]),
        'SDG': sorted([i for i in all_ids if i.startswith('SDG')]),
    }

    for group_name, ids in groups.items():
        if not ids:
            print(f"\n[{group_name}] no specimens, skipping")
            continue
        print(f"\n{'-'*50}")
        print(f"[{group_name}] {len(ids)} specimens")
        print(f"{'-'*50}")
        records = build_specimen_data(df, ids)
        if not records:
            print(f"[{group_name}] no valid data, skipping")
            continue
        html_str = generate_html(records, group_name, BANDWIDTH)
        out_path = os.path.join(OUT_DIR, f"kde_sphere_{group_name}.html")
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(html_str)
        size_mb = os.path.getsize(out_path) / (1024 * 1024)
        print(f"\n  OK  {out_path}")
        print(f"      {len(records)} specimens, {size_mb:.2f} MB")

    print(f"\n{'='*60}")
    print(f"Done. Output directory: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
