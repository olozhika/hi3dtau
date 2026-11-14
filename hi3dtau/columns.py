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


def background_fraction(n_los, grid, distance_kpc):
    """Fraction of the line-of-sight HI column lying behind ``distance_kpc``.

    Replaces the Li & Goldsmith (2003) smooth-disk erfc estimate and its
    empirical 0.7 floor with the column ratio of the 3D reconstruction.

    ``p_mean`` is the ratio of the mean columns, so that
    ``p_mean * N_total_mean == N_background_mean``; ``p_std`` is the scatter
    of the per-sample ratios.
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
        "p": p,
        "p_mean": float(p_mean),
        "p_std": float(p_std),
        "distance_kpc": float(distance_kpc),
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
