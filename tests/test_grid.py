import numpy as np
import healpy as hp
import pytest

from hi3dtau.grid import NSIDE, KPC_TO_CM, RadialGrid, ang2pix, unit_vector
from hi3dtau.columns import columns, p_3d, background_split_at


def test_ang2pix_matches_healpy():
    l = [0.0, 90.0, 180.0, 240.0, 359.9]
    b = [0.0, 0.0, 0.0, -15.77, 89.0]
    got = ang2pix(l, b)
    want = hp.ang2pix(NSIDE, np.radians(90 - np.array(b)), np.radians(l), nest=True)
    assert np.array_equal(got, want)


def test_unit_vector_orthonormal_directions():
    u = unit_vector([0.0, 90.0, 0.0], [0.0, 0.0, 90.0])
    assert np.allclose(u[0], [1, 0, 0])
    assert np.allclose(u[1], [0, 1, 0])
    assert np.allclose(u[2], [0, 0, 1])


def test_split_conserves_length():
    grid = RadialGrid([0.5, 1.5, 2.5], [0.0, 1.0, 2.0, 3.0])
    for d in [-1.0, 0.0, 0.25, 1.0, 1.7, 3.0, 10.0]:
        front, behind = grid.split(d)
        assert np.allclose(front + behind, grid.widths)


def test_partial_bin_fraction():
    grid = RadialGrid([0.5, 1.5], [0.0, 1.0, 2.0])
    front, behind = grid.split(0.5)
    assert np.allclose(front, [0.5, 0.0])
    assert np.allclose(behind, [0.5, 1.0])


def test_uniform_density_column():
    edges = np.array([0.0, 1.0, 2.0, 4.0])
    grid = RadialGrid(0.5 * (edges[1:] + edges[:-1]), edges)
    n = np.full((3, 2), 1.0)
    out = columns(n, grid, 1.5)
    assert np.isclose(out["N_total_mean"], 4.0 * KPC_TO_CM)
    assert np.isclose(out["N_foreground_mean"], 1.5 * KPC_TO_CM)
    assert np.isclose(out["N_background_mean"], 2.5 * KPC_TO_CM)


def test_p_3d_bounds():
    edges = np.array([0.0, 1.0, 2.0, 4.0])
    grid = RadialGrid(0.5 * (edges[1:] + edges[:-1]), edges)
    n = np.full((3, 2), 1.0)
    assert np.isclose(p_3d(n, grid, -1.0)["p_3d_mean"], 1.0)
    assert np.isclose(p_3d(n, grid, 10.0)["p_3d_mean"], 0.0)
    assert np.isclose(p_3d(n, grid, 1.5)["p_3d_mean"], 2.5 / 4.0)


def test_background_split_consistent():
    edges = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
    grid = RadialGrid(0.5 * (edges[1:] + edges[:-1]), edges)
    n = np.full((4, 2), 1.0)
    out = background_split_at(n, grid, 0.3, radius_kpc=1.25)
    total = out["N_background_within_radius"] + out["N_background_beyond_radius"]
    ref = columns(n, grid, 0.3)["N_background"]
    assert np.allclose(total, ref)
    assert np.allclose(out["N_background_beyond_radius"], (4.0 - 1.25) * KPC_TO_CM)
