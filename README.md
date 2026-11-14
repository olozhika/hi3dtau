# hi3dtau

HI foreground/background columns, the background fraction `p`, and HI optical
depths for any line of sight, read directly from the 3D reconstruction of the
Milky Way by [Söding et al. (2025)](https://doi.org/10.1051/0004-6361/202451361)
(A&A 693, A139).

For HINSA work this replaces the smooth-disk `erfc` estimate of Li & Goldsmith
(2003) with the column split of an actual 3D
map. Every quantity is reported as the mean over the authors' 8 posterior
samples together with the sample scatter.

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
p = hi3dtau.p(120.0, 30.0, 1.0)
print(p["p_mean"], p["p_std"])          # fraction of HI behind the cloud

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
    m.p(120.0, 30.0, 1.0)
    m.background_split(120.0, 30.0, 1.0, radius_kpc=1.25)
    m.velocity_los(120.0, 30.0)              # km/s, (n_radial, n_samples)
    m.spectrum(120.0, 30.0, [-10, 0, 10])    # K, (n_velocities, n_samples)
```

## API

| call | returns |
| --- | --- |
| `columns(l, b, D)` | `N_foreground`, `N_background`, `N_total` in cm^-2, each with per-sample array and `_mean` / `_std` |
| `p(l, b, D)` | `p` per sample, `p_mean`, `p_std`; `p_mean = N_background_mean / N_total_mean` |
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
| `samples_densities_hpixr.fits` | 1.2 GB | columns, `p` |
| `samples_auxiliary_hpixr.fits` | 3.1 GB | `tau`, `spectrum`, `velocity_los` |

Cache directory: `$HI3DTAU_DATA`, else a `data/` folder next to the package
if it already holds the files, else `~/.cache/hi3dtau`.

Nothing has to be downloaded by hand. To fetch everything up front (for a batch
job, or to check the download before a long run):

```python
import hi3dtau

hi3dtau.download()              # both products, MD5-checked
hi3dtau.download("densities")   # only what columns and p need
```

```bash
hi3dtau --download              # same, from the command line
```

Downloads resume after an interruption (a `.part` file is kept) and are
verified against the published MD5 on completion.

The native grid is Sun-centred HEALPix (Nside = 64, nested, ~55 arcmin) times
388 log-spaced radial bins from 0.05 to 30.2 kpc. Densities are in cm^-3 and
the files are memory-mapped, so only the accessed pixels are read.

