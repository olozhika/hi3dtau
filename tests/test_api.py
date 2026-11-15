"""End-to-end checks; skipped when the Zenodo data are not present."""

import os

import numpy as np
import pytest

from hi3dtau import data as _data
from hi3dtau.api import HI3D

pytestmark = pytest.mark.skipif(
    not (
        os.path.exists(_data.data_path("densities"))
        and os.path.getsize(_data.data_path("densities")) == _data.FILES["densities"][1]
    ),
    reason="density data not downloaded",
)


@pytest.fixture(scope="module")
def model():
    with HI3D() as m:
        yield m


def test_columns_positive_and_ordered(model):
    out = model.columns(120.0, 30.0, 1.2)
    assert out["N_total_mean"] > 0
    assert 0 <= out["N_foreground_mean"] < out["N_total_mean"]
    assert np.allclose(
        out["N_foreground_mean"] + out["N_background_mean"], out["N_total_mean"],
        rtol=1e-6,
    )


def test_p_3d_monotonic_in_distance(model):
    vals = [model.p_3d(120.0, 30.0, d)["p_3d_mean"] for d in (0.2, 1.0, 3.0, 8.0)]
    assert np.all(np.diff(vals) <= 1e-6)


def test_p_3d_matches_column_ratio(model):
    cols = model.columns(120.0, 30.0, 1.0)
    p = model.p_3d(120.0, 30.0, 1.0)["p_3d_mean"]
    assert np.isclose(p, cols["N_background_mean"] / cols["N_total_mean"], rtol=1e-6)


def test_p_v_positive_means_total(model):
    pv = model.p_v(120.0, 30.0, 1.0, velocity_kms=0.0)
    assert 0.0 <= pv["p_v_mean"] <= 1.0
    t = model.tau(120.0, 30.0, 1.0, 0.0)
    assert np.isclose(
        pv["p_v_mean"],
        t["tau_background_mean"] / (t["tau_foreground_mean"] + t["tau_background_mean"]),
        rtol=1e-6,
    )


def test_background_split_sums_to_background(model):
    split = model.background_split(120.0, 30.0, 1.0, radius_kpc=1.25)
    cols = model.columns(120.0, 30.0, 1.0)
    assert np.isclose(
        split["N_background_within_radius_mean"]
        + split["N_background_beyond_radius_mean"],
        cols["N_background_mean"],
        rtol=1e-6,
    )


def test_density_units_physical(model):
    n = model.density(120.0, 30.0)
    assert n.shape == (388, 8)
    assert np.nanmedian(n) > 0
    assert np.nanmax(n) < 1e4
