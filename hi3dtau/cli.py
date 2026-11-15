"""Command-line interface: one line of sight, one report."""

from __future__ import annotations

import argparse
import json
import sys

from .api import HI3D
from .data import download


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="hi3dtau",
        description=(
            "Foreground/background HI columns and optical depths from the "
            "Soeding et al. (2025) 3D Milky Way reconstruction."
        ),
    )
    ap.add_argument("--l", type=float, help="Galactic longitude [deg]")
    ap.add_argument("--b", type=float, help="Galactic latitude [deg]")
    ap.add_argument("--d", type=float, help="distance [kpc]")
    ap.add_argument(
        "--download",
        nargs="?",
        const="all",
        choices=["all", "densities", "auxiliary"],
        help="download the Zenodo data products and exit",
    )
    ap.add_argument("--v", type=float, default=None, help="LSR velocity [km/s] for tau")
    ap.add_argument("--spin", type=float, default=200.0, help="HI spin temperature [K]")
    ap.add_argument("--radius", type=float, default=1.25, help="reference radius [kpc]")
    ap.add_argument("--data-dir", default=None, help="directory with the FITS files")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    args = ap.parse_args(argv)

    if args.download is not None:
        download(None if args.download == "all" else args.download, args.data_dir)
        print("hi3dtau: data ready", file=sys.stderr)
        return 0

    if args.l is None or args.b is None or args.d is None:
        ap.error("--l, --b and --d are required unless --download is used")

    with HI3D(data_dir=args.data_dir) as model:
        cols = model.columns(args.l, args.b, args.d)
        frac = model.p_3d(args.l, args.b, args.d)
        split = model.background_split(args.l, args.b, args.d, args.radius)
        result = {
            "l_deg": args.l,
            "b_deg": args.b,
            "distance_kpc": args.d,
            "p_3d": frac["p_3d_mean"],
            "p_3d_std": frac["p_3d_std"],
            "N_foreground": cols["N_foreground_mean"],
            "N_foreground_std": cols["N_foreground_std"],
            "N_background": cols["N_background_mean"],
            "N_background_std": cols["N_background_std"],
            "N_total": cols["N_total_mean"],
            "N_total_std": cols["N_total_std"],
            f"fraction_background_beyond_{args.radius:g}kpc": split[
                "fraction_beyond_radius_mean"
            ],
            f"N_background_beyond_{args.radius:g}kpc": split[
                "N_background_beyond_radius_mean"
            ],
        }
        if args.v is not None:
            t = model.tau(args.l, args.b, args.d, args.v, T_spin=args.spin)
            pv = model.p_v(args.l, args.b, args.d, args.v, T_spin=args.spin)
            result.update(
                {
                    "velocity_kms": args.v,
                    "T_spin": args.spin,
                    "p_v": pv["p_v_mean"],
                    "p_v_std": pv["p_v_std"],
                    "tau_foreground": t["tau_foreground_mean"],
                    "tau_foreground_std": t["tau_foreground_std"],
                    "tau_background": t["tau_background_mean"],
                    "tau_background_std": t["tau_background_std"],
                }
            )

    if args.json:
        print(json.dumps(result, indent=2))
        return 0

    print(f"l = {args.l:.3f} deg, b = {args.b:.3f} deg, D = {args.d:.3f} kpc")
    print(f"p_3d (background fraction) = {result['p_3d']:.3f} +/- {result['p_3d_std']:.3f}")
    print(
        f"N_HI = {result['N_total']:.3e} +/- {result['N_total_std']:.1e} cm^-2 "
        f"(foreground {result['N_foreground']:.2e}, "
        f"background {result['N_background']:.2e})"
    )
    print(
        f"background column beyond {args.radius:g} kpc: "
        f"{100 * result[f'fraction_background_beyond_{args.radius:g}kpc']:.0f}%"
    )
    if args.v is not None:
        print(
            f"p_v(v={args.v:g} km/s) = {result['p_v']:.3f} +/- {result['p_v_std']:.3f}"
        )
        print(
            f"tau(v={args.v:g} km/s, T_spin={args.spin:g} K) = "
            f"{result['tau_foreground']:.3f} +/- {result['tau_foreground_std']:.3f} "
            f"(foreground), "
            f"{result['tau_background']:.3f} +/- {result['tau_background_std']:.3f} "
            f"(background)"
        )
    return 0


if __name__ == "__main__":
    main()
