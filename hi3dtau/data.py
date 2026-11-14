"""Access to the Soeding et al. (2025) 3D HI/CO reconstruction data products.

Data are published on Zenodo under DOI 10.5281/zenodo.12578443 and are
downloaded on first use. Two products are used:

``samples_densities_hpixr.fits``
    HI and H2 number density, 8 posterior samples.
``samples_auxiliary_hpixr.fits``
    Galactic velocity field (3 components) and HI/CO line widths,
    8 posterior samples.

Both live on a Sun-centred HEALPix (Nside=64, nested) x log(radius) grid.
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
import urllib.error
import urllib.request

ZENODO_RECORD = "12578443"
ZENODO_URL = "https://zenodo.org/records/{record}/files/{name}?download=1"

FILES = {
    "densities": (
        "samples_densities_hpixr.fits",
        1220564160,
        "e9012b22a1ff2df55100a61a02fec6f1",
    ),
    "auxiliary": (
        "samples_auxiliary_hpixr.fits",
        3051388800,
        "ca04d4a265c10f7297cce2f28d8dfe5b",
    ),
}

_NSIDE = 64


def _repo_data_dir() -> str:
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))


def default_data_dir() -> str:
    env = os.environ.get("HI3DTAU_DATA")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    repo = _repo_data_dir()
    if any(os.path.exists(os.path.join(repo, f[0])) for f in FILES.values()):
        return repo
    return os.path.join(
        os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache")), "hi3dtau"
    )


def data_path(kind: str, data_dir: str | None = None) -> str:
    if kind not in FILES:
        raise ValueError(f"unknown data product {kind!r}, choose from {sorted(FILES)}")
    return os.path.join(data_dir or default_data_dir(), FILES[kind][0])


def _md5(path: str, block: int = 8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(block)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _fmt_bytes(n: float) -> str:
    if n >= 1e9:
        return f"{n / 1e9:.2f} GB"
    if n >= 1e6:
        return f"{n / 1e6:.0f} MB"
    return f"{n / 1e3:.0f} kB"


def _progress(done: int, total: int, rate: float, t0: float) -> None:
    eta = (total - done) / rate if rate > 1e5 else float("nan")
    eta_txt = f"{int(eta // 60)}:{int(eta % 60):02d}" if eta == eta else "--:--"
    print(
        f"hi3dtau:   {_fmt_bytes(done)}/{_fmt_bytes(total)} "
        f"({100 * done / total:.0f}%, {rate / 1e6:.1f} MB/s, ETA {eta_txt})",
        file=sys.stderr,
    )


def _download(name: str, dest: str, expected_size: int, retries: int = 3) -> None:
    """Fetch one product into ``dest``, resuming a partial ``.part`` file."""
    url = ZENODO_URL.format(record=ZENODO_RECORD, name=name)
    part = dest + ".part"
    for attempt in range(1, retries + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have > expected_size:
            os.remove(part)
            have = 0
        if have == expected_size:
            os.replace(part, dest)
            return
        headers = {"Range": f"bytes={have}-"} if have else {}
        req = urllib.request.Request(url, headers=headers)
        if have:
            print(
                f"hi3dtau: resuming {name} at {_fmt_bytes(have)}", file=sys.stderr
            )
        else:
            print(f"hi3dtau: downloading {name} -> {dest}", file=sys.stderr)
        t0 = time.time()
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                status = getattr(resp, "status", None) or resp.getcode()
                if have and status != 206:
                    have = 0
                    if os.path.exists(part):
                        os.remove(part)
                    print(
                        "hi3dtau:   server ignored the range request, restarting",
                        file=sys.stderr,
                    )
                with open(part, "ab") as out:
                    done = have
                    next_report = (done // (200 * 1024 * 1024) + 1) * 200 * 1024 * 1024
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)
                        done += len(chunk)
                        if done >= next_report:
                            _progress(done, expected_size,
                                      (done - have) / max(time.time() - t0, 1e-6), t0)
                            next_report += 200 * 1024 * 1024
                    out.flush()
                    os.fsync(out.fileno())
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            if attempt == retries:
                raise
            wait = 2 ** attempt
            print(
                f"hi3dtau:   {exc}; retrying in {wait} s ({attempt}/{retries})",
                file=sys.stderr,
            )
            time.sleep(wait)
            continue
        got = os.path.getsize(part)
        if got != expected_size:
            if attempt == retries:
                raise IOError(
                    f"incomplete download for {name}: got {got} of {expected_size} bytes"
                )
            continue
        _progress(expected_size, expected_size,
                  (expected_size - have) / max(time.time() - t0, 1e-6), t0)
        os.replace(part, dest)
        return


def ensure(kind: str, data_dir: str | None = None, verify: bool = False) -> str:
    """Return the local path of a data product, downloading it if needed.

    A fresh download is always checked against the published MD5. An existing
    file is trusted on its size unless ``verify=True``, which re-hashes it.
    """
    name, size, md5 = FILES[kind]
    dest = data_path(kind, data_dir)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        if verify:
            print(f"hi3dtau: verifying {name}", file=sys.stderr)
            if _md5(dest) != md5:
                raise IOError(
                    f"checksum mismatch for {dest}; delete it to download again"
                )
        return dest
    _download(name, dest, size)
    print(f"hi3dtau: verifying {name}", file=sys.stderr)
    if _md5(dest) != md5:
        os.remove(dest)
        raise IOError(
            f"checksum mismatch for newly downloaded {dest}; "
            "the file was removed, please retry"
        )
    return dest


def download(kind: str | None = None, data_dir: str | None = None) -> dict:
    """Download one product (``"densities"`` or ``"auxiliary"``) or both.

    Useful to fetch everything before a long batch job. Returns a mapping from
    product name to local path.
    """
    kinds = [kind] if kind else list(FILES)
    out = {}
    for k in kinds:
        if k not in FILES:
            raise ValueError(f"unknown data product {k!r}, choose from {sorted(FILES)}")
        out[k] = ensure(k, data_dir)
    return out
