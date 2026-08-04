"""
perturb_spharm.py
=================
Synthetic annotation-error engine for the SP-SPHARM robustness analysis
(analysis/robustness/annotation_perturbation/).

Answers reviewer 2: manual placement of direction landmarks is itself error-prone
(Liu et al. 2026). This turns that limitation statement into quantitative evidence
by perturbing the UPSTREAM UNIT VECTORS and re-running the entire
vector -> vMF-sKDE -> DH grid -> spherical-harmonic -> power-spectrum chain, exactly
as the production pipeline does. Noise is NEVER added to the power spectra directly.

Two perturbations:
  P1 polarity : flip a fraction f of scar vectors (proximal/distal mislabelled)
  P2 angle    : isotropic random rotation of each unit vector, s.d. sigma degrees

A scar-dropout perturbation was tried and REMOVED. The sparsest EXP specimen has 10
scars, so even 20% dropout leaves 8 and the 3-scar floor never once triggered: the
perturbation only ever operated in the data-rich regime, i.e. it asked "what if a
specimen with plenty of scars loses a few?", whose answer is necessarily "nothing
much". The question that matters — how few scars can a specimen have — is a
downsampling analysis (k = 3/5/8/10/15) and belongs in its own script.

For every replicate the engine recomputes SPI and the fabric metrics (E, I) from the
SAME perturbed vector array (see Engine.run: U is perturbed once per specimen, then
the KDE, SPI and E/I are all derived from that one array), so the three methods are
compared on byte-identical data with no extra sampling noise.

Fabric is EXACTLY invariant under P1 — the orientation tensor sum(u u') is unchanged
by u -> -u (Mark 1973) — which makes the polarity contrast a mathematical identity
rather than an empirical finding. P2 is the real experiment: no method has analytic
immunity to small rotations, so the degradation ordering is genuinely unknown a priori.

REUSES the production implementation verbatim:
  SPHARM_modules.spherical_kde.make_sphere_grid / fit_vmf_kde   (h = 0.35, kappa = 8.16)
  SPHARM_modules.kde_to_spharm.compute_spharm_features           (lmax = 20)
  SPHARM_modules.power_spectrum.compute_power_spectrum
Parameters match kde_to_spharm_main.py: BANDWIDTH 0.35, 72x36 source grid,
DH_SIZE 64 (64 x 128), LMAX 20.

ONE OPTIMISATION, mathematically identity-preserving: kde_vector_to_dh_grid()
rebuilds a fixed (n_lat, n_lon, n_src) interpolation-weight tensor on every call.
That tensor depends only on the two fixed grids, so it is hoisted out of the loop
and pre-normalised once into a (n_lat*n_lon, n_src) linear operator. The per-
specimen result is the same arithmetic in a different association order; the
--verify mode checks the reconstructed baseline spectra against the committed
SPHARM_direction.csv to prove it.

Usage
-----
  python perturb_spharm.py --verify              # reproduce baseline, print max error
  python perturb_spharm.py --timing              # time one full 58-specimen pass
  python perturb_spharm.py --run --reps 100      # full sweep -> replicate CSVs
"""

import argparse
import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/project/analysis/scripts")

from SPHARM_modules.spherical_kde import make_sphere_grid, fit_vmf_kde
from SPHARM_modules.kde_to_spharm import compute_spharm_features

# ============================================================
# Config — identical to kde_to_spharm_main.py
# ============================================================
DERIVED_DIR = "/project/analysis/data/derived_data"
OUT_DIR     = "/project/analysis/robustness/annotation_perturbation"
SRC_CSV     = f"{DERIVED_DIR}/directions_aligned_svd.csv"
REF_CSV     = f"{DERIVED_DIR}/SPHARM_direction.csv"

BANDWIDTH = 0.35
N_BEARING = 72
N_PLUNGE  = 36
LMAX      = 20
DH_SIZE   = 64
SEED      = 42

# Perturbation levels
POLARITY_LEVELS = [0.02, 0.05, 0.10, 0.15, 0.20]
ANGLE_LEVELS    = [5.0, 10.0, 15.0, 20.0]

POWER_COLS_KEEP = [f"power_l{l}" for l in range(1, 7)]   # SP-SPHARM uses l = 1-6


# ============================================================
# Fixed interpolation operator (hoisted out of the per-specimen loop)
# ============================================================
def build_dh_operator(G: np.ndarray, dh_size: int = DH_SIZE):
    """
    Pre-normalised DH interpolation operator.

    Reproduces kde_to_spharm.kde_vector_to_dh_grid()'s weighting exactly:
        w = exp(50 * clip(dot(t, s), -1, 1)),  grid = (w @ kde) / rowsum(w)
    but builds w once, since it depends only on the source grid G and the DH grid.

    Returns (Wn, sin_colat) with Wn shape (n_lat * n_lon, n_src), rows summing to 1.
    """
    n_lat, n_lon = dh_size, 2 * dh_size
    colat_dh = np.linspace(0, np.pi,     n_lat, endpoint=False)
    lon_dh   = np.linspace(0, 2*np.pi,   n_lon, endpoint=False)
    TH, PH   = np.meshgrid(colat_dh, lon_dh, indexing="ij")

    T = np.column_stack([(np.sin(TH) * np.cos(PH)).ravel(),
                         (np.sin(TH) * np.sin(PH)).ravel(),
                         np.cos(TH).ravel()])                  # (n_lat*n_lon, 3)

    # Source grid unit vectors, in the same bearing/plunge convention the
    # production code derives from G.
    bearing   = np.arctan2(G[:, 1], G[:, 0])
    plunge    = np.arcsin(np.clip(G[:, 2], -1, 1))
    colat_src = np.pi / 2 - plunge
    S = np.column_stack([np.sin(colat_src) * np.cos(bearing),
                         np.sin(colat_src) * np.sin(bearing),
                         np.cos(colat_src)])                   # (n_src, 3)

    dot = np.clip(T @ S.T, -1, 1)
    W   = np.exp(50.0 * dot)
    Wn  = W / W.sum(axis=1, keepdims=True)
    return Wn, np.sin(colat_dh)


def kde_to_power(kde_mat: np.ndarray, Wn: np.ndarray, sin_colat: np.ndarray,
                 lmax: int = LMAX) -> np.ndarray:
    """
    KDE vectors -> DH grids -> SH expansion -> normalised power spectra.
    kde_mat: (n_spec, n_src).  Returns (n_spec, lmax+1) normalised power.
    """
    n_lat = len(sin_colat)
    n_lon = Wn.shape[0] // n_lat
    grids = kde_mat @ Wn.T                                     # (n_spec, n_lat*n_lon)
    out = np.empty((kde_mat.shape[0], lmax + 1))
    for i in range(kde_mat.shape[0]):
        g = np.clip(grids[i].reshape(n_lat, n_lon), 0, None)
        area = (g * sin_colat[:, None]).sum()
        g = g / (area if area > 0 else 1.0)
        out[i] = compute_spharm_features(g, lmax=lmax)["norm_power"]
    return out


# ============================================================
# Direction metrics — verbatim ports of spharm_analysis.R:214-231
# ============================================================
def compute_spi(U: np.ndarray) -> float:
    """spharm_analysis.R:214-216 — resultant length / n (polarity SENSITIVE)."""
    return float(np.sqrt((U.sum(axis=0) ** 2).sum()) / len(U))


def compute_ei(U: np.ndarray):
    """spharm_analysis.R:218-231 — orientation tensor (polarity INSENSITIVE)."""
    T = (U.T @ U) / len(U)
    lam = np.sort(np.linalg.eigvalsh(T))[::-1]
    lam = np.maximum(lam, 0.0)
    if lam[0] <= 1e-10:
        return np.nan, np.nan
    return float(1.0 - lam[1] / lam[0]), float(lam[2] / lam[0])


# ============================================================
# Perturbations — applied to the UPSTREAM UNIT VECTORS
# ============================================================
def perturb_polarity(U: np.ndarray, f: float, rng) -> np.ndarray:
    """P1: flip a random fraction f of vectors (u -> -u)."""
    if f <= 0:
        return U
    n = len(U)
    k = int(round(f * n))
    if k == 0:
        return U
    idx = rng.choice(n, size=k, replace=False)
    V = U.copy()
    V[idx] *= -1.0
    return V


def perturb_angle(U: np.ndarray, sigma_deg: float, rng) -> np.ndarray:
    """
    P2: isotropic random rotation of each unit vector.

    Implementation: for each vector u, draw a uniformly random direction in u's
    TANGENT PLANE (via Gram-Schmidt on a random 3-vector) and rotate u by an angle
    theta about the axis perpendicular to both. theta ~ |N(0, sigma)|, with the
    tangent direction uniform on [0, 2*pi), so the perturbation is isotropic about
    u and the result stays exactly on the unit sphere (Rodrigues rotation).
    """
    if sigma_deg <= 0:
        return U
    n = len(U)
    sigma = np.deg2rad(sigma_deg)

    # Random tangent direction: project a random Gaussian onto u's tangent plane.
    R = rng.normal(size=(n, 3))
    R -= (np.sum(R * U, axis=1, keepdims=True)) * U
    nrm = np.linalg.norm(R, axis=1, keepdims=True)
    bad = (nrm[:, 0] < 1e-12)
    if bad.any():                       # degenerate draw: retry those rows once
        R2 = rng.normal(size=(bad.sum(), 3))
        R2 -= (np.sum(R2 * U[bad], axis=1, keepdims=True)) * U[bad]
        R[bad] = R2
        nrm = np.linalg.norm(R, axis=1, keepdims=True)
    Tdir = R / np.maximum(nrm, 1e-300)

    theta = rng.normal(loc=0.0, scale=sigma, size=(n, 1))
    # Rotate u toward Tdir by theta; u and Tdir are orthonormal, so this is exact.
    V = np.cos(theta) * U + np.sin(theta) * Tdir
    return V / np.linalg.norm(V, axis=1, keepdims=True)


# ============================================================
# One full pass: perturbed vectors -> spectra + SPI + fabric
# ============================================================
class Engine:
    def __init__(self):
        self.G = make_sphere_grid(N_BEARING, N_PLUNGE)
        self.Wn, self.sin_colat = build_dh_operator(self.G)
        df = pd.read_csv(SRC_CSV)
        df = df[df["ID"].str.startswith("EXP")].dropna(subset=["ux", "uy", "uz"])
        self.ids = sorted(df["ID"].unique())
        self.groups = {i: df[df["ID"] == i][["ux", "uy", "uz"]].to_numpy(float)
                       for i in self.ids}
        self.typ = {i: df[df["ID"] == i]["Typology"].iloc[0] for i in self.ids}
        self.kappa = 1.0 / BANDWIDTH ** 2

    def _kde(self, U: np.ndarray) -> np.ndarray:
        # identical to fit_vmf_kde(); inlined so the fixed grid is reused
        dot = self.G @ U.T
        dens = np.mean(np.exp(self.kappa * dot), axis=1)
        return dens / dens.sum()

    def run(self, perturb=None):
        """perturb: callable(ID, U) -> U' (or None for the baseline)."""
        kde = np.empty((len(self.ids), len(self.G)))
        spi = np.empty(len(self.ids))
        Evals = np.empty(len(self.ids))
        Ivals = np.empty(len(self.ids))
        floored = 0
        for j, i in enumerate(self.ids):
            U = self.groups[i]
            if perturb is not None:
                U, fl = perturb(i, U)
                floored += int(fl)
            kde[j] = self._kde(U)
            spi[j] = compute_spi(U)
            Evals[j], Ivals[j] = compute_ei(U)
        power = kde_to_power(kde, self.Wn, self.sin_colat)
        out = pd.DataFrame({"ID": self.ids,
                            "Typology": [self.typ[i] for i in self.ids],
                            "SPI": spi, "E": Evals, "I": Ivals})
        for l in range(LMAX + 1):
            out[f"power_l{l}"] = power[:, l]
        return out, floored


# ============================================================
# Modes
# ============================================================
def do_verify(eng: Engine):
    base, _ = eng.run(None)
    ref = pd.read_csv(REF_CSV)
    ref = ref[ref["ID"].isin(base["ID"])].set_index("ID").loc[base["ID"]]
    cols = [f"power_l{l}" for l in range(1, 21)]
    diff = np.abs(base.set_index("ID")[cols].to_numpy() - ref[cols].to_numpy())
    print(f"specimens compared : {len(base)}")
    print(f"max |power diff|   : {diff.max():.3e}")
    print(f"max over l = 1-6   : {diff[:, :6].max():.3e}")
    ok = diff[:, :6].max() < 1e-6
    print("VERIFY:", "PASS — fast operator reproduces the committed spectra"
          if ok else "FAIL — do not trust the perturbation results")
    return ok


def do_timing(eng: Engine):
    t0 = time.time(); eng.run(None); t1 = time.time()
    rng = np.random.default_rng(0)
    t2 = time.time()
    eng.run(lambda i, U: (perturb_polarity(U, 0.10, rng), False))
    t3 = time.time()
    print(f"baseline pass (58 specimens) : {t1-t0:6.2f} s")
    print(f"perturbed pass               : {t3-t2:6.2f} s")
    n_cond = len(POLARITY_LEVELS) + len(ANGLE_LEVELS)
    for reps in (50, 100):
        print(f"  python total, R={reps:3d}, {n_cond} conditions : "
              f"{(t3-t2)*n_cond*reps/60:6.1f} min")


def do_run(eng: Engine, reps: int):
    os.makedirs(OUT_DIR, exist_ok=True)
    base, _ = eng.run(None)
    base.insert(0, "rep", 0); base.insert(0, "level", 0.0)
    base.insert(0, "perturbation", "baseline")
    frames = [base]

    def seed_for(kind, level, rep):
        # crc32 is stable across processes; Python's hash() on str is salted by
        # PYTHONHASHSEED and would make the sweep irreproducible between runs.
        key = f"{kind}|{float(level):.4f}|{int(rep)}".encode()
        return (zlib.crc32(key) ^ SEED) & 0x7FFFFFFF

    t0 = time.time()
    for kind, levels in (("polarity", POLARITY_LEVELS),
                         ("angle",    ANGLE_LEVELS)):
        for lv in levels:
            for r in range(1, reps + 1):
                rng = np.random.default_rng(seed_for(kind, lv, r))
                if kind == "polarity":
                    fn = lambda i, U, rng=rng, lv=lv: (perturb_polarity(U, lv, rng), False)
                else:
                    fn = lambda i, U, rng=rng, lv=lv: (perturb_angle(U, lv, rng), False)
                out, _ = eng.run(fn)
                out.insert(0, "rep", r); out.insert(0, "level", lv)
                out.insert(0, "perturbation", kind)
                frames.append(out)
            print(f"  {kind} {lv}: {reps} reps done "
                  f"({time.time()-t0:.0f}s elapsed)", flush=True)

    allout = pd.concat(frames, ignore_index=True)
    keep = ["perturbation", "level", "rep", "ID", "Typology", "SPI", "E", "I"] + POWER_COLS_KEEP
    allout[keep].to_csv(f"{OUT_DIR}/perturbed_descriptors.csv", index=False)
    print(f"\nWrote perturbed_descriptors.csv ({len(allout)} rows)")
    print(f"Total python time: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--timing", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--reps", type=int, default=100)
    a = ap.parse_args()
    eng = Engine()
    print(f"EXP specimens: {len(eng.ids)}; total scars: "
          f"{sum(len(v) for v in eng.groups.values())}")
    if a.verify:
        sys.exit(0 if do_verify(eng) else 1)
    if a.timing:
        do_timing(eng)
    if a.run:
        do_run(eng, a.reps)
