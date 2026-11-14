"""Geometry of the Sun-centred HEALPix x log(radius) grid."""

from __future__ import annotations

import numpy as np
import healpy as hp

NSIDE = 64
NPIX = 12 * NSIDE * NSIDE
KPC_TO_CM = 3.0856775814913673e21


def ang2pix(l_deg, b_deg):
    """HEALPix (nested) pixel index for Galactic longitude/latitude in degrees."""
    l = np.asarray(l_deg, dtype=float)
    b = np.asarray(b_deg, dtype=float)
    return hp.ang2pix(NSIDE, np.radians(90.0 - b), np.radians(l), nest=True)


def unit_vector(l_deg, b_deg):
    """Unit vector in the Sun-centred Galactic frame (x to GC, z to NGP)."""
    l = np.radians(np.asarray(l_deg, dtype=float))
    b = np.radians(np.asarray(b_deg, dtype=float))
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1)


class RadialGrid:
    """Log-spaced radial bins, with exact partial-bin splitting at a distance."""

    def __init__(self, centres_kpc, edges_kpc):
        self.centres = np.asarray(centres_kpc, dtype=float)
        self.edges = np.asarray(edges_kpc, dtype=float)
        if self.centres.size != self.edges.size - 1:
            raise ValueError("centres and edges are inconsistent")
        self.widths = np.diff(self.edges)
        self.r_min = float(self.edges[0])
        self.r_max = float(self.edges[-1])

    def split(self, distance_kpc):
        """Path lengths (kpc) of every bin in front of and behind a distance."""
        d = float(distance_kpc)
        lo = self.edges[:-1]
        hi = self.edges[1:]
        front = np.clip(np.minimum(hi, d) - lo, 0.0, None)
        behind = np.clip(hi - np.maximum(lo, d), 0.0, None)
        return front, behind


def los_density(image, ipix):
    """Density along one line of sight, shaped (n_radial, n_samples)."""
    return np.asarray(image[ipix], dtype=float)
