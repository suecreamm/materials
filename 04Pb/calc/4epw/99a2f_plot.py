#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EPW alpha^2F plotter.

Features
--------
- Auto-detects *.a2f.* files in '.' then '..'
- Safely reads the numeric alpha^2F block before the EPW summary section
- Reads phonon-smearing values and integrated lambda from the summary
- Saves PNG + PDF
- Can plot one smearing curve or all smearing curves
- Styled to resemble a dark-background publication/presentation figure

Examples
--------
  python a2f_plot.py
  python a2f_plot.py --file pb.a2f.01.0.075
  python a2f_plot.py --smearing 0.10
  python a2f_plot.py --all-smearings
  python a2f_plot.py --unit mev
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SEARCH_DIRS = [Path("."), Path("..")]
CM1_PER_MEV = 8.06554429
CM1_PER_THZ = 33.3564095


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Plot EPW Eliashberg spectral function alpha^2F."
    )
    ap.add_argument(
        "--file",
        default=None,
        help="EPW *.a2f.* file. If omitted, auto-detected in '.' then '..'."
    )
    ap.add_argument(
        "--smearing",
        type=float,
        default=None,
        help="Phonon smearing in meV to plot. Default: first available value."
    )
    ap.add_argument(
        "--all-smearings",
        action="store_true",
        help="Plot every phonon-smearing curve."
    )
    ap.add_argument(
        "--unit",
        choices=["mev", "thz", "cm-1"],
        default="mev",
        help="Frequency-axis unit (default: meV)."
    )
    ap.add_argument(
        "--out",
        default="a2f",
        help="Output basename (default: a2f)."
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution (default: 300)."
    )
    ap.add_argument(
        "--title",
        default="",
        help="Optional figure title."
    )
    ap.add_argument("--xmin", type=float, default=None)
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("--ymax", type=float, default=None)
    return ap.parse_args()


def resolve_a2f(file_arg: Optional[str]) -> Path:
    if file_arg:
        p = Path(file_arg)
        if p.is_file():
            return p
        raise FileNotFoundError(f"File not found: {file_arg}")

    patterns = ("*.a2f.*", "*.a2f")
    for d in SEARCH_DIRS:
        for pattern in patterns:
            hits = sorted(
                p for p in d.glob(pattern)
                if p.is_file() and "a2f_tr" not in p.name
            )
            if hits:
                return hits[0]

    raise FileNotFoundError(
        "No EPW alpha^2F file found in '.' or '..'. "
        "Use --file pb.a2f.01.0.075"
    )


def _numbers_from_line(line: str) -> List[float]:
    vals = []
    for tok in line.replace("#", " ").split():
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    return vals


def load_epw_a2f(
    path: Path
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[float]]:
    """
    Returns
    -------
    omega_mev : (n,)
        Frequency grid from EPW. For these EPW files the first column is in meV.
    a2f : (n, nsmear)
        alpha^2F curves for each phonon smearing.
    smearings_mev : (nsmear,)
        Phonon-smearing values from EPW summary.
    lambdas : (nsmear,)
        Integrated electron-phonon coupling for each smearing.
    lambda_sum : float or None
        EPW 'Summed el-ph coupling' value.
    """
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    numeric_rows: List[List[float]] = []
    summary_start = len(lines)

    for i, line in enumerate(lines):
        if "Integrated el-ph coupling" in line:
            summary_start = i
            break

        s = line.strip()
        if not s or s.startswith("#"):
            continue

        vals = _numbers_from_line(s)
        if len(vals) >= 2:
            if numeric_rows and len(vals) != len(numeric_rows[0]):
                break
            numeric_rows.append(vals)

    if not numeric_rows:
        raise ValueError(f"No alpha^2F numeric block found in {path}")

    ncol = len(numeric_rows[0])
    arr = np.asarray([row for row in numeric_rows if len(row) == ncol], dtype=float)

    omega_mev = arr[:, 0]
    a2f = arr[:, 1:]

    smearings: List[float] = []
    lambdas: List[float] = []
    lambda_sum: Optional[float] = None

    for i in range(summary_start, len(lines)):
        line = lines[i]

        if "Integrated el-ph coupling" in line and i + 1 < len(lines):
            lambdas = _numbers_from_line(lines[i + 1])

        if "Phonon smearing" in line and i + 1 < len(lines):
            smearings = _numbers_from_line(lines[i + 1])

        if "Summed el-ph coupling" in line:
            vals = _numbers_from_line(line)
            if vals:
                lambda_sum = vals[-1]

    ncurve = a2f.shape[1]

    if not smearings:
        smearings = list(range(1, ncurve + 1))
        print("[WARN] Phonon-smearing values not found; using curve indices.")

    if not lambdas:
        lambdas = [np.nan] * ncurve
        print("[WARN] Integrated lambda values not found.")

    smearings_arr = np.asarray(smearings[:ncurve], dtype=float)
    lambdas_arr = np.asarray(lambdas[:ncurve], dtype=float)

    if len(smearings_arr) != ncurve:
        raise ValueError(
            f"alpha^2F has {ncurve} curves but only {len(smearings_arr)} "
            "smearing values were found."
        )

    return omega_mev, a2f, smearings_arr, lambdas_arr, lambda_sum


def convert_x(x_mev: np.ndarray, unit: str) -> Tuple[np.ndarray, str]:
    if unit == "mev":
        return x_mev, "Phonon energy (meV)"
    if unit == "thz":
        return x_mev * CM1_PER_MEV / CM1_PER_THZ, "Frequency (THz)"
    if unit == "cm-1":
        return x_mev * CM1_PER_MEV, r"Frequency (cm$^{-1}$)"
    raise ValueError(unit)


def choose_curve(smearings: np.ndarray, requested: Optional[float]) -> int:
    if requested is None:
        return 0

    idx = int(np.argmin(np.abs(smearings - requested)))
    if not np.isclose(smearings[idx], requested, atol=1e-8):
        print(
            f"[WARN] Requested smearing {requested:g} meV not found. "
            f"Using nearest value {smearings[idx]:g} meV."
        )
    return idx


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("black")

    # ticks
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=16,
        colors="white",
        width=1.5,
        length=7,
        direction="out",
        top=True,
        right=True
    )

    # spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("white")
        spine.set_linewidth(1.2)

    # grid
    ax.grid(
        True,
        which="major",
        color="0.7",
        linewidth=1.8,
        alpha=1.0
    )


def make_single_curve_legend(ax: plt.Axes, sm: float, lam: float) -> None:
    label = rf"$\sigma_{{ph}}$ = {sm:g} meV"
    if np.isfinite(lam):
        label += "\n" + rf"$\lambda$ = {lam:.6f}"

    handles = ax.plot([], [], lw=2.5, label=label)
    leg = ax.legend(
        handles=handles,
        loc="upper left",
        fontsize=20,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        borderpad=0.4,
        handlelength=2.6,
        handletextpad=0.9,
        labelspacing=0.3
    )

    frame = leg.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_linewidth(1.5)

    for text in leg.get_texts():
        text.set_color("black")


def make_all_smearings_legend(ax: plt.Axes) -> None:
    leg = ax.legend(
        loc="upper left",
        fontsize=13,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        borderpad=0.4,
        handlelength=2.3,
        handletextpad=0.7,
        labelspacing=0.25
    )

    frame = leg.get_frame()
    frame.set_facecolor("white")
    frame.set_edgecolor("black")
    frame.set_linewidth(1.3)

    for text in leg.get_texts():
        text.set_color("black")


def main() -> None:
    args = parse_args()
    path = resolve_a2f(args.file)

    omega_mev, a2f, smearings, lambdas, lambda_sum = load_epw_a2f(path)
    x, xlabel = convert_x(omega_mev, args.unit)

    print(f"[INFO] alpha2F file : {path}")
    print(f"[INFO] points       : {len(x)}")
    print(f"[INFO] smearings    : {', '.join(f'{v:g}' for v in smearings)} meV")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    fig.patch.set_facecolor("black")
    style_axes(ax)

    if args.all_smearings:
        for j, sm in enumerate(smearings):
            lam = lambdas[j] if j < len(lambdas) else np.nan

            label = rf"$\sigma_{{ph}}$ = {sm:g} meV"
            if np.isfinite(lam):
                label += rf", $\lambda$ = {lam:.6f}"

            ax.plot(
                x,
                a2f[:, j],
                lw=2.0,
                label=label
            )

        make_all_smearings_legend(ax)
        print(f"[INFO] plotting     : all {len(smearings)} smearing curves")

    else:
        j = choose_curve(smearings, args.smearing)
        sm = smearings[j]
        lam = lambdas[j] if j < len(lambdas) else np.nan

        # actual curve
        ax.plot(x, a2f[:, j], lw=2.4)

        # filled area
        ax.fill_between(x, 0.0, a2f[:, j], alpha=0.90)

        # separate boxed legend like the example image
        make_single_curve_legend(ax, sm, lam)

        print(f"[INFO] smearing     : {sm:g} meV")
        if np.isfinite(lam):
            print(f"[INFO] lambda       : {lam:.7f}")

    if lambda_sum is not None:
        print(f"[INFO] summed lambda: {lambda_sum:.7f}")

    # labels
    ax.set_xlabel(xlabel, fontsize=20, color="white", labelpad=10)
    ax.set_ylabel(r"$\alpha^2F(\omega)$", fontsize=20, color="white", labelpad=10)

    # optional title only
    if args.title:
        ax.set_title(args.title, fontsize=20, color="white", pad=12)

    # zero line
    ax.axhline(0.0, lw=1.0, color="white")

    # limits
    if args.xmin is not None or args.xmax is not None:
        ax.set_xlim(
            args.xmin if args.xmin is not None else float(np.nanmin(x)),
            args.xmax if args.xmax is not None else float(np.nanmax(x)),
        )
    else:
        ax.set_xlim(float(np.nanmin(x)), float(np.nanmax(x)))

    ax.set_ylim(bottom=0.0)
    if args.ymax is not None:
        ax.set_ylim(0.0, args.ymax)

    fig.tight_layout(pad=1.0)

    out_base = (
        str(Path(args.out).with_suffix(""))
        if Path(args.out).suffix.lower() in (".png", ".pdf")
        else args.out
    )
    out_png = f"{out_base}.png"
    out_pdf = f"{out_base}.pdf"

    fig.savefig(
        out_png,
        dpi=args.dpi,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight"
    )
    fig.savefig(
        out_pdf,
        facecolor=fig.get_facecolor(),
        bbox_inches="tight"
    )
    plt.close(fig)

    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")


if __name__ == "__main__":
    main()
