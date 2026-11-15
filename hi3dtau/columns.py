"""Foreground/background HI columns along a line of sight."""

from __future__ import annotations

import numpy as np

from .grid import KPC_TO_CM


def _summarise(values, axis=0):
    mean = np.mean(values, axis=axis)
    std = np.std(values, axis=axis, ddof=1) if values.shape[axis] > 1 else np.zeros_like(mean)
    return mean, std


def columns(n_los, grid, distance_kpc):
    """HI columns in cm^-2 in front of and behind ``distance_kpc``.

    Parameters
    ----------
    n_los : ndarray (n_radial, n_samples)
        HI number density in cm^-3 along one line of sight.
    grid : RadialGrid
    distance_kpc : float

    Returns a dict with per-sample arrays and their mean/std over samples.
    """
    front, behind = grid.split(distance_kpc)
    n_fg = (n_los * front[:, None]).sum(axis=0) * KPC_TO_CM
    n_bg = (n_los * behind[:, None]).sum(axis=0) * KPC_TO_CM
    n_tot = n_fg + n_bg
    fg_mean, fg_std = _summarise(n_fg)
    bg_mean, bg_std = _summarise(n_bg)
    tot_mean, tot_std = _summarise(n_tot)
    return {
        "N_foreground": n_fg,
        "N_background": n_bg,
        "N_total": n_tot,
        "N_foreground_mean": float(fg_mean),
        "N_foreground_std": float(fg_std),
        "N_background_mean": float(bg_mean),
        "N_background_std": float(bg_std),
        "N_total_mean": float(tot_mean),
        "N_total_std": float(tot_std),
        "distance_kpc": float(distance_kpc),
        "r_min_kpc": grid.r_min,
        "r_max_kpc": grid.r_max,
    }


def p_3d(n_los, grid, distance_kpc):
    """Fraction of the total line-of-sight HI column behind ``distance_kpc``.

    Velocity-independent: every radial bin contributes regardless of its
    line-of-sight velocity. Replaces the Li & Goldsmith (2003) smooth-disk
    erfc estimate and its empirical 0.7 floor with the column ratio of the
    3D reconstruction.

    ``p_3d_mean`` is the ratio of the mean columns, so that
    ``p_3d_mean * N_total_mean == N_background_mean``; ``p_3d_std`` is the
    scatter of the per-sample ratios.
    """
    front, behind = grid.split(distance_kpc)
    n_fg = (n_los * front[:, None]).sum(axis=0)
    n_bg = (n_los * behind[:, None]).sum(axis=0)
    n_tot = n_fg + n_bg
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(n_tot > 0, n_bg / n_tot, np.nan)
    _, p_std = _summarise(p)
    p_mean = np.mean(n_bg) / np.mean(n_tot)
    return {
        "p_3d": p,
        "p_3d_mean": float(p_mean),
        "p_3d_std": float(p_std),
        "distance_kpc": float(distance_kpc),
    }


def p_v(n_los, sigma_los, v_los, grid, distance_kpc, velocity_kms, T_spin=200.0):
    """Fraction of the 21 cm optical depth at ``velocity_kms`` behind the cloud.

    Velocity-resolved counterpart of :func:`p_3d`: each radial bin is weighted
    by its Gaussian line profile overlap with ``velocity_kms``, then split
    front/back at ``distance_kpc``. In the optically thin limit this is the
    fraction of the line emission at that velocity available as a background.

    ``p_v_mean`` is the ratio of the mean foreground/background optical depths;
    ``p_v_std`` is the scatter of the per-sample ratios. ``tau_*`` are the
    per-sample optical depths.
    """
    from .tau import tau as _tau
    from .tau import summarise as _summarise_tau

    front, behind = grid.split(distance_kpc)
    tau_fg = _tau(n_los, sigma_los, v_los, front, velocity_kms, T_spin=T_spin)
    tau_bg = _tau(n_los, sigma_los, v_los, behind, velocity_kms, T_spin=T_spin)
    tau_tot = tau_fg + tau_bg
    fg_mean, fg_std = _summarise_tau(tau_fg)
    bg_mean, bg_std = _summarise_tau(tau_bg)
    tot_mean, tot_std = _summarise_tau(tau_tot)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(tau_tot > 0, tau_bg / tau_tot, np.nan)
    _, p_std = _summarise(p)
    p_mean = np.mean(tau_bg) / np.mean(tau_tot)
    return {
        "p_v": p,
        "p_v_mean": float(p_mean),
        "p_v_std": float(p_std),
        "tau_foreground": tau_fg,
        "tau_background": tau_bg,
        "tau_total": tau_tot,
        "tau_foreground_mean": float(fg_mean),
        "tau_foreground_std": float(fg_std),
        "tau_background_mean": float(bg_mean),
        "tau_background_std": float(bg_std),
        "tau_total_mean": float(tot_mean),
        "tau_total_std": float(tot_std),
        "distance_kpc": float(distance_kpc),
        "velocity_kms": float(velocity_kms),
        "T_spin": float(T_spin),
    }


def background_split_at(n_los, grid, distance_kpc, radius_kpc=1.25):
    """Background column inside/outside a reference radius (default 1.25 kpc).

    The 1.25 kpc scale is where Soeding et al. (2025) validate their
    reconstruction against 3D dust maps. This reports how much of the
    background column rests on the less-validated far field.
    """
    lo = grid.edges[:-1]
    hi = grid.edges[1:]
    d = float(distance_kpc)
    near = np.clip(np.minimum(hi, radius_kpc) - np.maximum(lo, d), 0.0, None)
    far = np.clip(hi - np.maximum(np.maximum(lo, d), radius_kpc), 0.0, None)
    n_near = (n_los * near[:, None]).sum(axis=0) * KPC_TO_CM
    n_far = (n_los * far[:, None]).sum(axis=0) * KPC_TO_CM
    n_bg = n_near + n_far
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_far = np.where(n_bg > 0, n_far / n_bg, np.nan)
    frac_mean, frac_std = _summarise(frac_far)
    near_mean, near_std = _summarise(n_near)
    far_mean, far_std = _summarise(n_far)
    return {
        "N_background_within_radius": n_near,
        "N_background_beyond_radius": n_far,
        "N_background_within_radius_mean": float(near_mean),
        "N_background_within_radius_std": float(near_std),
        "N_background_beyond_radius_mean": float(far_mean),
        "N_background_beyond_radius_std": float(far_std),
        "fraction_beyond_radius": frac_far,
        "fraction_beyond_radius_mean": float(frac_mean),
        "fraction_beyond_radius_std": float(frac_std),
        "radius_kpc": float(radius_kpc),
    }
