"""Minimal standalone example.

Runs without installing the package: the package sits next to this script''s
parent directory, so the directory is put on ``sys.path`` when needed. The Zenodo data products (~4.3 GB total) are downloaded on first use
into ``~/.cache/hi3dtau`` (or the directory given by ``HI3DTAU_DATA``).
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if (_ROOT / "hi3dtau" / "__init__.py").exists():
    sys.path.insert(0, str(_ROOT))

import hi3dtau  # noqa: E402

L, B, D = 120.0, 30.0, 1.0  # deg, deg, kpc
V = 0.0                     # km/s, line-centre velocity

res = hi3dtau.p_3d(L, B, D)
print(f"p_3d = {res['p_3d_mean']:.3f} +/- {res['p_3d_std']:.3f}")

pv = hi3dtau.p_v(L, B, D, velocity_kms=V, T_spin=200.0)
print(f"p_v   = {pv['p_v_mean']:.3f} +/- {pv['p_v_std']:.3f}  (at v={V:g} km/s)")

cols = hi3dtau.columns(L, B, D)
print(f"N_fg = {cols['N_foreground_mean']:.3e} cm^-2")
print(f"N_bg = {cols['N_background_mean']:.3e} cm^-2")
print(f"N_HI = {cols['N_total_mean']:.3e} cm^-2")

t = hi3dtau.tau(L, B, D, velocity_kms=V, T_spin=200.0)
print(f"tau_fg = {t['tau_foreground_mean']:.3f}")
print(f"tau_bg = {t['tau_background_mean']:.3f}")
