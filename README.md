# hi3dtau

HI foreground/background columns, the background fractions `p_3d` / `p_v`, and
HI optical depths for any line of sight, read directly from the 3D
reconstruction of the Milky Way by
[Söding et al. (2025)](https://doi.org/10.1051/0004-6361/202451361)
(A&A 693, A139).

For HINSA work this replaces the smooth-disk `erfc` estimate of Li & Goldsmith
(2003) with the split of an actual 3D map. Every quantity is reported as the
mean over the authors' 8 posterior samples together with the sample scatter.

## `p_3d` vs `p_v`

Two background fractions, computed the same way (column/optical depth behind
the cloud divided by the total) but differing in what counts as "behind":

| name | definition | notes |
| --- | --- | --- |
| `p_3d(l, b, D)` | **all** of the line-of-sight HI column split at distance `D`, ignoring velocity: `p_3d = N_bg / (N_fg + N_bg)` | velocity-agnostic; robust, cheap; the "3D" generalisation of Li & Goldsmith's `erfc` |
| `p_v(l, b, D, velocity_kms)` | same spatial split, but only the HI whose line profile overlaps `velocity_kms`: `p_v = tau_bg(v) / (tau_fg(v) + tau_bg(v))` | the physically relevant quantity for HINSA (self-absorption happens in one velocity channel); needs the auxiliary velocity field, so it shares the caveats of `tau` |

In the optically thin limit the optical depth is proportional to the line
emission, so `p_v` is the fraction of the 21 cm emission at the cloud's
velocity that comes from behind it. The two differ when much of the background
column lies at velocities far from the cloud.

## Install

Copy the `hi3dtau` folder to wherever you want to use it, and do `import hi3dtau` in python

Only Python >= 3.9 and `numpy`, `astropy`, `healpy` are needed.

Optional, to use it from any directory:

```bash
pip install .
# or, for development
pip install -e ".[dev]"
```

## Quickstart

```python
import hi3dtau

# Galactic coordinates in degrees, distance in kpc
p3d = hi3dtau.p_3d(120.0, 30.0, 1.0)
print(p3d["p_3d_mean"], p3d["p_3d_std"])   # fraction of total HI behind the cloud

# Same split, restricted to one velocity channel (physically right for HINSA)
pv = hi3dtau.p_v(120.0, 30.0, 1.0, velocity_kms=0.0)
print(pv["p_v_mean"], pv["p_v_std"])       # fraction of the 21 cm emission at
                                           # v=0 km/s coming from behind

cols = hi3dtau.columns(120.0, 30.0, 1.0)
print(cols["N_foreground_mean"], cols["N_background_mean"], cols["N_total_mean"])

# Optical depth at the line-centre velocity (50-100 K is typical for cold HI)
t = hi3dtau.tau(120.0, 30.0, 1.0, velocity_kms=0.0, T_spin=100.0)
print(t["tau_foreground_mean"], t["tau_background_mean"])
```

`python examples/quickstart.py` runs the same demo as a script.

Reuse a loaded map instead of reopening the FITS files on every call:

```python
from hi3dtau import HI3D

with HI3D() as m:
    m.columns(120.0, 30.0, 1.0)
    m.p_3d(120.0, 30.0, 1.0)
    m.p_v(120.0, 30.0, 1.0, velocity_kms=0.0)
    m.background_split(120.0, 30.0, 1.0, radius_kpc=1.25)
    m.velocity_los(120.0, 30.0)              # km/s, (n_radial, n_samples)
    m.spectrum(120.0, 30.0, [-10, 0, 10])    # K, (n_velocities, n_samples)
```

## API

| call | returns |
| --- | --- |
| `columns(l, b, D)` | `N_foreground`, `N_background`, `N_total` in cm^-2, each with per-sample array and `_mean` / `_std` |
| `p_3d(l, b, D)` | `p_3d`, `p_3d_mean`, `p_3d_std`; fraction of the total HI column behind `D`, velocity-agnostic |
| `p_v(l, b, D, velocity_kms, T_spin=200)` | `p_v`, `p_v_mean`, `p_v_std` plus `tau_foreground/background/total`; fraction of the 21 cm opacity at `velocity_kms` behind `D` |
| `background_split(l, b, D, radius_kpc=1.25)` | background column inside/outside a reference radius and the fraction beyond it |
| `velocity_los(l, b)` | reconstructed line-of-sight velocity [km/s] |
| `tau(l, b, D, velocity_kms, T_spin=200)` | foreground and background HI optical depth at that velocity |
| `spectrum(l, b, velocities_kms, T_spin=200)` | model HI brightness temperature spectrum [K] |
| `summary(...)` | all of the above as one dictionary |

Command line:

```bash
hi3dtau --l 120 --b 30 --d 1.0 --v 0 --spin 100
hi3dtau --l 120 --b 30 --d 1.0 --json
```

## Data

The two posterior-sample FITS files are downloaded from
[Zenodo 10.5281/zenodo.12578443](https://doi.org/10.5281/zenodo.12578443) on
first use and verified against their MD5 checksums:

| file | size | loaded when |
| --- | --- | --- |
| `samples_densities_hpixr.fits` | 1.2 GB | columns, `p_3d` |
| `samples_auxiliary_hpixr.fits` | 3.1 GB | `tau`, `spectrum`, `velocity_los` |

Cache directory: `$HI3DTAU_DATA`, else a `data/` folder next to the package
if it already holds the files, else `~/.cache/hi3dtau`.

Nothing has to be downloaded by hand. To fetch everything up front (for a batch
job, or to check the download before a long run):

```python
import hi3dtau

hi3dtau.download()              # both products, MD5-checked
hi3dtau.download("densities")   # only what columns and p_3d need
```

```bash
hi3dtau --download              # same, from the command line
```

Downloads resume after an interruption (a `.part` file is kept) and are
verified against the published MD5 on completion.

The native grid is Sun-centred HEALPix (Nside = 64, nested, ~55 arcmin) times
388 log-spaced radial bins from 0.05 to 30.2 kpc. Densities are in cm^-3 and
the files are memory-mapped, so only the accessed pixels are read.

## Caveats

- **Spin temperature.** The reconstruction fixes `T_spin = 200 K`. Optical
  depth scales as `1 / T_spin`; pass a lower `T_spin` for cold foreground HI.
- **Resolution.** ~55 arcmin angularly; radial bins run from ~1 pc near the Sun
  to ~130 pc at 8 kpc. Columns and `p_3d` are line-of-sight integrals and are the
  robust products, but cloud-scale structure is not resolved.
- **Velocity field.** The published `VELOCITY IN X/Y/Z-DIRECTION` maps are not
  in the standard Galactic Cartesian basis: the Y axis is mirrored. This
  package therefore projects them as `v_los = u_x v_x - u_y v_y + u_z v_z` with
  `u = (cos b cos l, cos b sin l, sin b)`. The maps are relative to the Sun, and
  this projection is the one that reproduces observed CO velocities, whereas the
  naive all-plus projection does not.
- **Velocity-resolved optical depth is fragile.** At 55 arcmin the reconstructed
  velocity field is a beam average. When `|v_los(D) - v_0|` exceeds ~10 km/s the
  derived `tau` becomes very small and is no longer meaningful. Compare
  `spectrum` with the observed profile before using `tau`. Columns and `p_3d` do
  not depend on this.

## Validation

For 29 high-latitude sightlines (`|b| > 10 deg`), the total HI column from the
map agrees with the optically thin HI4PI column: median
`N_3D / N_HI4PI = 1.04`, 16-84% range 0.95-1.10. At lower latitudes the
thin-limit estimate is expected to fall below the true column.

## Tests

From the downloaded folder:

```bash
pytest -q
```

## Citation

If you use this package, cite the data paper and the data set:

- Söding, L., Edenhofer, G., Enßlin, T. A., et al. 2025, A&A 693, A139
- Zenodo data: doi:10.5281/zenodo.12578443

## Troubleshooting

`AttributeError: module 'hi3dtau' has no attribute 'p'`

Python imported the *downloaded folder* `hi3dtau/` (an empty namespace package)
instead of the package inside it. This happens when Python is started one level
**above** the downloaded folder instead of inside it:

```bash
cd hi3dtau            # the folder containing pyproject.toml
python
```

To check which copy was imported:

```python
import hi3dtau
print(hi3dtau.__file__)   # None means the folder was imported, not the package
```

Alternatively `pip install .` and use it from anywhere.

## License

MIT. See `LICENSE`.
