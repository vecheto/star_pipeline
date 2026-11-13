"""
Relative photometry with an ensemble comparison star (method from ../brown_dwarfs/metodo.md).

Arrays are (n_frames, n_stars); magnitudes are instrumental.
"""
import numpy as np
from astropy.stats import sigma_clip

POGSON = 2.5 / np.log(10)


def to_mag(flux, flux_err=None):
    """m = -2.5 log10(F),  sigma_m = 1.0857 sigma_F / F.  Non-positive flux -> NaN."""
    flux = np.where(flux > 0, flux, np.nan)
    mag = -2.5 * np.log10(flux)
    mag_err = None if flux_err is None else POGSON * flux_err / flux
    return mag, mag_err


def complete_stars(mag):
    """Boolean mask of stars measured in every frame (brown_dwarfs: clean_stars)."""
    return np.isfinite(mag).all(axis=0)


def ensemble_star(mag, mag_err=None, sigma=3):
    """
    Ensemble star per frame from the reference stars (target excluded by the caller).

        delta_i(t) = m_i(t) - median_t(m_i)
    Per frame, delta_i is rejected if |delta_i - mean(delta)| > sigma * std(delta), iterating
    until no rejections (point by point: a star can be rejected in one frame and used in another).
    With N <= 10 references nothing can be rejected at 3 sigma: max deviation is sqrt(N-1) std.
        mag_ens = sum(w_i delta_i) / sum(w_i),  w_i = 1 / sigma_i(t)^2
        mag_ens_err = (sum w_i)^(-1/2)
    Without errors (isophotal) the mean is unweighted and mag_ens_err is NaN.

    Returns mag_ens, mag_ens_err, rejected (bool, same shape as mag).
    """
    delta = mag - np.nanmedian(mag, axis=0)
    clipped = sigma_clip(delta, sigma=sigma, maxiters=None, cenfunc='mean', stdfunc='std', axis=1)
    keep = ~np.ma.getmaskarray(clipped) & np.isfinite(delta)

    weights = np.ones_like(delta) if mag_err is None else 1 / mag_err**2
    weights = np.where(keep, weights, 0)
    sum_w = weights.sum(axis=1)
    mag_ens = np.nansum(weights * np.where(keep, delta, 0), axis=1) / sum_w
    mag_ens_err = np.full_like(mag_ens, np.nan) if mag_err is None else 1 / np.sqrt(sum_w)
    return mag_ens, mag_ens_err, ~keep & np.isfinite(delta)


def relative_phot(mag_t, mag_t_err, mag_ens, mag_ens_err):
    """mag_relative = m_target - mag_ens;  error = sqrt(sigma_target^2 + sigma_ens^2)."""
    rel = mag_t - mag_ens
    rel_err = None if mag_t_err is None else np.sqrt(mag_t_err**2 + mag_ens_err**2)
    return rel, rel_err
