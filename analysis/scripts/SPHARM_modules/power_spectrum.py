# power_spectrum.py
# Per-degree power spectrum of a spherical-harmonic expansion.
# Used by SPHARM_main.py (morphology) and kde_to_spharm.py (direction).

import numpy as np


def compute_power_spectrum(clm_sh, lmax: int = None) -> dict:
    """
    Per-degree power spectrum of a pyshtools SHCoeffs object.

    Parameters
    ----------
    clm_sh : pyshtools.SHCoeffs
    lmax   : int, optional. If None, uses all available degrees.

    Returns
    -------
    dict with keys:
        raw_power  : ndarray (lmax+1,)  raw power per degree
        norm_power : ndarray (lmax+1,)  AC-normalised power (l=0 set to 0)
    """
    raw_power = clm_sh.spectrum()
    if lmax is not None:
        raw_power = raw_power[:lmax + 1]

    power_no_dc    = raw_power.copy()
    power_no_dc[0] = 0.0
    total_ac       = power_no_dc.sum()
    norm_power     = power_no_dc / total_ac if total_ac > 0 else power_no_dc

    return {
        "raw_power"  : raw_power,
        "norm_power" : norm_power,
    }
