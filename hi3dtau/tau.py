"""HI optical depth from the reconstructed number density and velocity field.

The absorption coefficient follows Soeding et al. (2025, Eq. 8):

    kappa(v) = 679.168 kpc^-1 (K / T_spin) (n_HI / cm^-3)
               (km/s / sigma_v) exp(-(v - v_los)^2 / 2 sigma_v^2)

with a Maxwellian line profile and the upper level assumed populated,
n_u/n_l = 3.
"""

from __future__ import annotations

import numpy as np

KAPPA_COEFF = 679.168


def kappa(n, sigma, v_los, velocity, T_spin=200.0):
    sigma = np.maximum(np.asarray(sigma, dtype=float), 1e-6)
    v_los = np.asarray(v_los, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    return (
        KAPPA_COEFF
        * (1.0 / float(T_spin))
        * np.asarray(n, dtype=float)
        / sigma
        * np.exp(-0.5 * ((velocity - v_los) / sigma) ** 2)
    )


def tau(n_los, sigma_los, v_los, lengths_kpc, velocity, T_spin=200.0):
    """Optical depth at ``velocity`` [km/s] for one line of sight.

    Parameters
    ----------
    n_los, sigma_los, v_los : ndarray (n_radial, n_samples)
        HI density [cm^-3], line width [km/s] and line-of-sight velocity
        [km/s] in each radial bin, per posterior sample.
    lengths_kpc : ndarray (n_radial,)
        Path length of each radial bin, or the part of it in front of / behind
        the cloud. Summing in front gives the foreground optical depth.
    velocity : float or array
        LSR velocity at which to evaluate the optical depth [km/s].
    """
    k = kappa(n_los, sigma_los, v_los, velocity, T_spin=T_spin)
    return (k * np.asarray(lengths_kpc, dtype=float)[:, None]).sum(axis=0)


def summarise(values):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    std = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    return mean, std


def brightness_temperature(n_los, sigma_los, v_los, widths_kpc, velocities, T_spin=200.0):
    """Model HI brightness temperature spectrum of one line of sight.

    Uses the discretised radiative transfer of Soeding et al. (2025, Eq. 14)
    with the source function fixed to ``T_spin``. Returns an array shaped
    (n_velocities, n_samples) in kelvin.
    """
    n = np.asarray(n_los, dtype=float)[:, None, :]
    sigma = np.asarray(sigma_los, dtype=float)[:, None, :]
    vl = np.asarray(v_los, dtype=float)[:, None, :]
    v = np.asarray(velocities, dtype=float)[None, :, None]
    k = kappa(n, sigma, vl, v, T_spin=T_spin)
    tau_bin = k * np.asarray(widths_kpc, dtype=float)[:, None, None]
    tau_front = np.cumsum(tau_bin, axis=0) - tau_bin
    emission = float(T_spin) * (1.0 - np.exp(-tau_bin)) * np.exp(-tau_front)
    return emission.sum(axis=0)

