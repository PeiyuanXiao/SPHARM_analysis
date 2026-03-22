import numpy as np
from scipy.stats import entropy as scipy_entropy


def compute_spectral_entropy(
    clm_sh,
    lmax: int = None,
) -> dict:
    """
    Compute SHE and spectral entropy from a pyshtools SHCoeffs object.

    Parameters
    ----------
    clm_sh : pyshtools.SHCoeffs
    lmax   : int, optional. If None, uses all available degrees.

    Returns
    -------
    dict with keys:
        raw_power        : ndarray (lmax+1,)  raw power per degree
        norm_power       : ndarray (lmax+1,)  normalised power (l=0 set to 0)
        SHE              : float              total power (sum of raw_power)
        spectral_entropy : float              Shannon entropy of norm_power[1:]
    """
    raw_power = clm_sh.spectrum()
    if lmax is not None:
        raw_power = raw_power[:lmax + 1]

    SHE = float(np.sum(raw_power))

    power_no_dc    = raw_power.copy()
    power_no_dc[0] = 0.0
    total_ac       = power_no_dc.sum()
    norm_power     = power_no_dc / total_ac if total_ac > 0 else power_no_dc

    p = norm_power[1:]
    p = p[p > 0]
    spectral_entropy = float(scipy_entropy(p))

    return {
        "raw_power"        : raw_power,
        "norm_power"       : norm_power,
        "SHE"              : SHE,
        "spectral_entropy" : spectral_entropy,
    }