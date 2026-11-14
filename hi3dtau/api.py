"""High-level interface to the 3D HI reconstruction."""

from __future__ import annotations

import os

import numpy as np
from astropy.io import fits

from . import data as _data
from .columns import (
    background_fraction,
    background_split_at,
    columns as _columns,
)
from .grid import RadialGrid, ang2pix, los_density, unit_vector
from .tau import brightness_temperature, summarise, tau as _tau

_IMAGE = (fits.ImageHDU, fits.PrimaryHDU)


def _image_hdus(hdul):
    return [h for h in hdul if isinstance(h, _IMAGE) and h.data is not None]


def _find(hdus, *needles, exclude=()):
    for h in hdus:
        name = str(h.name or h.header.get("EXTNAME", "")).upper()
        if all(n in name for n in needles) and not any(e in name for e in exclude):
            return h
    return None


class HI3D:
    """Foreground/background columns and optical depths at (l, b, D).

    Parameters
    ----------
    data_dir : str, optional
        Directory holding the Zenodo FITS files. Defaults to ``$HI3DTAU_DATA``,
        then ``hi3dtau/data`` in this repository, then ``~/.cache/hi3dtau``.
    """

    def __init__(self, data_dir=None, load_auxiliary=False):
        self.data_dir = data_dir
        dens_path = _data.ensure("densities", data_dir)
        self._dens = fits.open(dens_path, memmap=True)
        imgs = _image_hdus(self._dens)
        self._hi = _find(imgs, "HI", "DENS")
        self._h2 = _find(imgs, "H2") or _find(imgs, "H2", "DENS")
        if self._hi is None:
            raise RuntimeError(
                "HI density HDU not found, available: "
                + ", ".join(str(h.name) for h in imgs)
            )
        self._centres = np.asarray(self._dens[1].data.field(0), dtype=float)
        self._edges = np.asarray(self._dens[2].data.field(0), dtype=float)
        self.grid = RadialGrid(self._centres, self._edges)
        self._aux = None
        if load_auxiliary:
            self._load_auxiliary()

    def close(self):
        if self._aux is not None:
            self._aux.close()
            self._aux = None
        if self._dens is not None:
            self._dens.close()
            self._dens = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _load_auxiliary(self):
        if self._aux is not None:
            return
        path = _data.ensure("auxiliary", self.data_dir)
        self._aux = fits.open(path, memmap=True)
        imgs = _image_hdus(self._aux)
        self._vx = _find(imgs, "X-DIRECTION") or _find(imgs, "VELOCITY", "X")
        self._vy = _find(imgs, "Y-DIRECTION") or _find(imgs, "VELOCITY", "Y")
        self._vz = _find(imgs, "Z-DIRECTION") or _find(imgs, "VELOCITY", "Z")
        self._sig_hi = _find(imgs, "WIDTH", "HI") or _find(imgs, "SIGMA", "HI")
        self._sig_co = _find(imgs, "WIDTH", "CO") or _find(imgs, "SIGMA", "CO")
        missing = [
            n
            for n, h in [
                ("velocity x", self._vx),
                ("velocity y", self._vy),
                ("velocity z", self._vz),
                ("HI line width", self._sig_hi),
            ]
            if h is None
        ]
        if missing:
            raise RuntimeError(
                "could not identify "
                + ", ".join(missing)
                + "; available image HDUs: "
                + ", ".join(str(h.name) for h in imgs)
            )

    @property
    def n_samples(self):
        return int(self._hi.data.shape[-1])

    def _pixel(self, l_deg, b_deg):
        ipix = ang2pix(l_deg, b_deg)
        if np.ndim(ipix) == 0:
            return int(ipix)
        return ipix

    def density(self, l_deg, b_deg):
        """HI number density along the line of sight, (n_radial, n_samples)."""
        return los_density(self._hi.data, self._pixel(l_deg, b_deg))

    def columns(self, l_deg, b_deg, distance_kpc):
        return _columns(self.density(l_deg, b_deg), self.grid, distance_kpc)

    def p(self, l_deg, b_deg, distance_kpc):
        """Fraction of the HI column behind the cloud (3D equivalent of p)."""
        return background_fraction(
            self.density(l_deg, b_deg), self.grid, distance_kpc
        )

    def background_split(self, l_deg, b_deg, distance_kpc, radius_kpc=1.25):
        return background_split_at(
            self.density(l_deg, b_deg), self.grid, distance_kpc, radius_kpc
        )

    def _vlos(self, ipix, l_deg, b_deg):
        u = unit_vector(l_deg, b_deg)
        vx = np.asarray(self._vx.data[ipix], dtype=float)
        vy = np.asarray(self._vy.data[ipix], dtype=float)
        vz = np.asarray(self._vz.data[ipix], dtype=float)
        return u[0] * vx - u[1] * vy + u[2] * vz

    def velocity_los(self, l_deg, b_deg):
        """Reconstructed line-of-sight velocity field, (n_radial, n_samples)."""
        self._load_auxiliary()
        return self._vlos(self._pixel(l_deg, b_deg), l_deg, b_deg)

    def tau(self, l_deg, b_deg, distance_kpc, velocity_kms, T_spin=200.0):
        """HI optical depth in front of and behind the cloud at ``velocity_kms``.

        Returns a dict with foreground/background optical depths per posterior
        sample together with their mean and standard deviation.
        """
        self._load_auxiliary()
        ipix = self._pixel(l_deg, b_deg)
        n = self.density(l_deg, b_deg)
        v_los = self._vlos(ipix, l_deg, b_deg)
        sigma = np.asarray(self._sig_hi.data[ipix], dtype=float)
        front, behind = self.grid.split(distance_kpc)
        tau_fg = _tau(n, sigma, v_los, front, velocity_kms, T_spin=T_spin)
        tau_bg = _tau(n, sigma, v_los, behind, velocity_kms, T_spin=T_spin)
        fg_mean, fg_std = summarise(tau_fg)
        bg_mean, bg_std = summarise(tau_bg)
        return {
            "tau_foreground": tau_fg,
            "tau_background": tau_bg,
            "tau_foreground_mean": fg_mean,
            "tau_foreground_std": fg_std,
            "tau_background_mean": bg_mean,
            "tau_background_std": bg_std,
            "tau_total_mean": fg_mean + bg_mean,
            "velocity_kms": float(velocity_kms),
            "T_spin": float(T_spin),
            "distance_kpc": float(distance_kpc),
        }

    def spectrum(self, l_deg, b_deg, velocities_kms, T_spin=200.0):
        """Model HI brightness temperature spectrum, (n_velocities, n_samples).

        Useful to check the reconstruction against an observed HI spectrum
        before trusting the derived optical depths.
        """
        self._load_auxiliary()
        ipix = self._pixel(l_deg, b_deg)
        n = self.density(l_deg, b_deg)
        v_los = self._vlos(ipix, l_deg, b_deg)
        sigma = np.asarray(self._sig_hi.data[ipix], dtype=float)
        return brightness_temperature(
            n, sigma, v_los, self.grid.widths, velocities_kms, T_spin=T_spin
        )

    def summary(self, l_deg, b_deg, distance_kpc, velocity_kms=None, T_spin=200.0):
        out = {
            "l_deg": float(l_deg),
            "b_deg": float(b_deg),
        }
        cols = self.columns(l_deg, b_deg, distance_kpc)
        out.update({k: v for k, v in cols.items() if not np.ndim(v)})
        out["p"] = self.p(l_deg, b_deg, distance_kpc)["p_mean"]
        out["p_std"] = self.p(l_deg, b_deg, distance_kpc)["p_std"]
        split = self.background_split(l_deg, b_deg, distance_kpc)
        out.update({k: v for k, v in split.items() if not np.ndim(v)})
        if velocity_kms is not None:
            t = self.tau(l_deg, b_deg, distance_kpc, velocity_kms, T_spin=T_spin)
            out.update({k: v for k, v in t.items() if not np.ndim(v)})
        return out


_default = {}


def get(data_dir=None) -> HI3D:
    key = os.path.abspath(data_dir) if data_dir else None
    if key not in _default:
        _default[key] = HI3D(data_dir=data_dir)
    return _default[key]


def p(l_deg, b_deg, distance_kpc, data_dir=None):
    return get(data_dir).p(l_deg, b_deg, distance_kpc)


def columns(l_deg, b_deg, distance_kpc, data_dir=None):
    return get(data_dir).columns(l_deg, b_deg, distance_kpc)


def tau(l_deg, b_deg, distance_kpc, velocity_kms, T_spin=200.0, data_dir=None):
    return get(data_dir).tau(l_deg, b_deg, distance_kpc, velocity_kms, T_spin=T_spin)


def summary(l_deg, b_deg, distance_kpc, velocity_kms=None, T_spin=200.0, data_dir=None):
    return get(data_dir).summary(
        l_deg, b_deg, distance_kpc, velocity_kms=velocity_kms, T_spin=T_spin
    )
