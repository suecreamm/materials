#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI-only phonon dispersion + phonon DOS plotter for Quantum ESPRESSO
--------------------------------------------------------------------
- Headless (Matplotlib 'Agg'), saves PNG + PDF by default.
- Dispersion supports *.freq.gp (gnuplot table) and raw *.freq (&plot header).
- Auto high-symmetry (HS) labeling:
  * If --qpath is not given and ./qpath.in exists, it is used by default.
  * If qpath file is missing/unreadable -> skip HS labeling (no crash).
- Adds phonon DOS panel on the right (width ratio 5:1) when DOS exists.
  * If DOS file is missing/unreadable -> skip DOS panel (dispersion only, no crash).
- Accepts either:
  * a concrete file path via --freq, OR
  * a PREFIX via --freq (auto-detects common QE filenames).

Typical QE workflow:
  ph.x -> q2r.x -> matdyn.x (dispersion + DOS)

Example:
  python plot_phonon_dispersion.py --freq TiS2_monolayer
"""

from __future__ import annotations

import os
import re
import argparse
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

# ---- unit conversions ----
CM1_PER_THZ = 33.356
CM1_PER_MEV = 8.066


# ==========================================================
# Argument parsing
# ==========================================================
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Headless phonon dispersion + DOS plotter for QE (.freq.gp/.freq + .phdos)."
    )
    ap.add_argument(
        "--freq", required=True,
        help=(
            "Path to .freq.gp (preferred) or raw .freq file. "
            "If a file does not exist, it is treated as PREFIX and common filenames are auto-detected."
        )
    )
    ap.add_argument(
        "--dos", default=None,
        help=(
            "Optional phonon DOS file (2 columns: freq(cm^-1) DOS). "
            "If omitted, will auto-detect from PREFIX (if --freq is PREFIX). "
            "If not found, DOS is skipped."
        )
    )
    ap.add_argument(
        "--qpath", default=None,
        help=(
            "Optional q-path file (band form) for HS labels. "
            "Default: use ./qpath.in if present. If not found, HS labels are skipped."
        )
    )
    ap.add_argument(
        "--out", default="phonon_dispersion",
        help=(
            "Output basename (default: phonon_dispersion). "
            "If you pass a name with extension, extension is ignored; both .png and .pdf are produced."
        )
    )
    ap.add_argument("--emin", type=float, default=None, help="Y-axis min.")
    ap.add_argument("--emax", type=float, default=None, help="Y-axis max.")
    ap.add_argument("--dpi", type=int, default=300, help="PNG DPI.")
    ap.add_argument("--title", default="", help="Figure title.")
    ap.add_argument(
        "--unit", choices=["mev", "thz", "cm-1"], default="mev",
        help="Output unit for frequencies (default: meV)."
    )
    return ap.parse_args()


# ==========================================================
# File auto-detection
# ==========================================================
def _first_existing(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p.exists() and p.is_file():
            return p
    return None


def resolve_inputs(
    freq_arg: str, dos_arg: Optional[str], qpath_arg: Optional[str]
) -> Tuple[Path, Optional[Path], Optional[Path], str]:
    """
    Resolve input files.
    Missing qpath/dos are treated as None (skip later), not fatal.

    Returns
    -------
    freq_path : Path
        The resolved dispersion file (.freq.gp or .freq or table).
    dos_path : Optional[Path]
        The resolved DOS file (2-column) or None.
    qpath_path : Optional[Path]
        The resolved qpath file for HS labels or None.
    prefix : str
        Prefix inferred if freq_arg was treated as PREFIX; otherwise "".
    """
    freq_candidate = Path(freq_arg)
    prefix = ""

    # 1) Resolve dispersion file (fatal if not found)
    if freq_candidate.exists() and freq_candidate.is_file():
        freq_path = freq_candidate
    else:
        prefix = freq_arg
        guesses = [
            Path(f"{prefix}_phband.freq.gp"),
            Path(f"{prefix}_phband.freq"),
            Path(f"{prefix}.freq.gp"),
            Path(f"{prefix}.freq"),
            Path(f"{prefix}_dispersion.freq.gp"),
            Path(f"{prefix}_dispersion.freq"),
            Path(f"{prefix}_phband.freq.gp.dat"),
            Path(f"{prefix}_phband.freq.dat"),
        ]
        found = _first_existing(guesses)
        if found is None:
            raise FileNotFoundError(
                f"Could not find dispersion file. Tried: {[str(p) for p in guesses]}"
            )
        freq_path = found

    # 2) Resolve DOS file (non-fatal)
    dos_path: Optional[Path] = None
    if dos_arg:
        p = Path(dos_arg)
        if p.exists() and p.is_file():
            dos_path = p
        else:
            print(f"[WARN] --dos given but not found. Skipping DOS: {dos_arg}")
            dos_path = None
    else:
        if prefix:
            dos_guesses = [
                Path(f"{prefix}_phdos"),
                Path(f"{prefix}.phdos"),
                Path(f"{prefix}_phdos.dat"),
                Path(f"{prefix}.phdos.dat"),
            ]
            dos_path = _first_existing(dos_guesses)

    # 3) Resolve qpath file (non-fatal)
    qpath_path: Optional[Path] = None
    if qpath_arg:
        p = Path(qpath_arg)
        if p.exists() and p.is_file():
            qpath_path = p
        else:
            print(f"[WARN] --qpath given but not found. Skipping HS labels: {qpath_arg}")
            qpath_path = None
    else:
        default_qpath = Path("qpath.in")
        if default_qpath.exists() and default_qpath.is_file():
            qpath_path = default_qpath

    return freq_path, dos_path, qpath_path, prefix


def sanitize_out_basename(out_arg: str) -> str:
    """
    Convert --out to a basename (strip extension if user gave one).
    """
    p = Path(out_arg)
    if p.suffix.lower() in [".png", ".pdf"]:
        return str(p.with_suffix(""))
    return str(p)


# ==========================================================
# Helpers (parsers)
# ==========================================================
def _is_numeric_line(s: str) -> bool:
    s = s.strip()
    if not s or s.startswith("#"):
        return False
    return bool(re.match(r"^[\s]*[+\-]?\d", s))


def load_freq_gp(path_gp: Path) -> Tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(str(path_gp))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad .freq.gp table: {path_gp}")
    return data[:, 0], data[:, 1:]


def load_freq_raw(path_freq: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Parse matdyn.x raw .freq file (&plot header)."""
    with open(path_freq, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]
    m = re.match(r"^\s*&plot\s+nbnd=\s*([0-9]+)\s*,\s*nks=\s*([0-9]+)\s*/", lines[0])
    if not m:
        raise ValueError(f"Not a raw .freq header: {lines[0]}")
    nbnd, nks = int(m.group(1)), int(m.group(2))

    q_list, bands = [], []
    i = 1
    while i < len(lines) and len(q_list) < nks:
        toks = lines[i].split()
        if len(toks) != 3:
            raise ValueError(f"Unexpected q-vector line at line {i+1}")
        q_list.append(toks)
        i += 1
        cur: List[float] = []
        while i < len(lines) and len(cur) < nbnd and _is_numeric_line(lines[i]):
            cur.extend(np.fromstring(lines[i], sep=" ").tolist())
            i += 1
        bands.append(cur[:nbnd])

    q_idx = np.arange(len(q_list), dtype=float)
    return q_idx, np.array(bands, dtype=float)


def load_freq_table(path_any: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load generic whitespace table with columns:
      q  f1  f2  ...
    """
    data = np.loadtxt(str(path_any))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad frequency table: {path_any}")
    return data[:, 0], data[:, 1:]


def load_dos_2col(path_dos: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load DOS file as 2 columns:
      freq(cm^-1)  dos
    Any extra columns are ignored (first 2 used).
    """
    data = np.loadtxt(str(path_dos))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad DOS table (need >=2 cols): {path_dos}")
    f_cm1 = data[:, 0].astype(float)
    dos = data[:, 1].astype(float)
    idx = np.argsort(f_cm1)
    return f_cm1[idx], dos[idx]


def parse_qpath_labels(qpath_path: Path, n_points_data: int) -> Tuple[List[str], List[int]]:
    """
    Read HS labels from QE-style q-path input (band form) and infer segment tick indices.

    If parsing fails for any reason, returns empty lists (skip HS labels).
    """
    try:
        with open(qpath_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
        n_hsp = int(lines[0].split()[0])
        labels: List[str] = []
        for i in range(1, 1 + n_hsp):
            if i >= len(lines):
                break
            toks = lines[i].split()
            lbl = toks[-1] if len(toks) >= 4 else f"P{i}"
            if re.fullmatch(r"(?i)(g|gamma|Γ)", lbl):
                lbl = "Γ"
            labels.append(lbl)

        n_segments = max(0, len(labels) - 1)
        if n_segments == 0:
            idx = list(range(len(labels)))
        else:
            step = (n_points_data - 1) / n_segments
            idx = [int(round(i * step)) for i in range(n_segments + 1)]
        return labels, idx
    except Exception as e:
        print(f"[WARN] Failed to parse qpath labels from {qpath_path}. Skipping HS labels. ({e})")
        return [], []


def convert_units_from_cm1(y_cm1: np.ndarray, unit: str) -> Tuple[np.ndarray, str]:
    """
    Convert a frequency axis given in cm^-1 into requested unit.
    """
    if unit == "cm-1":
        return y_cm1, "Frequency (cm$^{-1}$)"
    if unit == "thz":
        return y_cm1 / CM1_PER_THZ, "Frequency (THz)"
    if unit == "mev":
        return y_cm1 / CM1_PER_MEV, "Energy (meV)"
    raise ValueError(f"Unknown unit: {unit}")


# ==========================================================
# Plot
# ==========================================================
def plot_dispersion_only(
    q_path: np.ndarray,
    freqs_cm1: np.ndarray,
    out_base: str,
    emin: Optional[float],
    emax: Optional[float],
    full_title: str,
    unit: str,
    labels: List[str],
    tick_idx: List[int],
    dpi: int = 300,
) -> None:
    """
    Dispersion-only plot (used when DOS is missing).
    Saves PNG + PDF.
    """
    y_disp, ylabel = convert_units_from_cm1(freqs_cm1, unit)

    fig, ax = plt.subplots(figsize=(8, 6))

    for j in range(y_disp.shape[1]):
        ax.plot(q_path, y_disp[:, j], "-", lw=2.0, color="tab:blue")
    ax.axhline(0.0, lw=0.8, color="black")

    if (emin is not None) or (emax is not None):
        lo = emin if emin is not None else float(np.nanmin(y_disp))
        hi = emax if emax is not None else float(np.nanmax(y_disp))
        ax.set_ylim(lo, hi)
    ax.set_xlim(float(np.min(q_path)), float(np.max(q_path)))

    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(full_title, fontsize=18)
    ax.tick_params(axis="y", labelsize=12)

    if labels and tick_idx and len(labels) == len(tick_idx):
        tick_pos = [q_path[i] for i in tick_idx]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(labels, fontsize=18)
        for x in tick_pos:
            ax.axvline(x, linestyle="--", linewidth=0.7, color="gray", alpha=0.6)
    else:
        ax.set_xticks([])

    ax.grid(alpha=0.3)
    fig.tight_layout()

    out_png = f"{out_base}.png"
    out_pdf = f"{out_base}.pdf"
    fig.savefig(out_png, dpi=dpi, transparent=True)
    fig.savefig(out_pdf, transparent=True)
    plt.close(fig)

    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")


def plot_dispersion_with_dos(
    q_path: np.ndarray,
    freqs_cm1: np.ndarray,
    dos_freq_cm1: np.ndarray,
    dos_val: np.ndarray,
    out_base: str,
    emin: Optional[float],
    emax: Optional[float],
    full_title: str,
    unit: str,
    labels: List[str],
    tick_idx: List[int],
    dpi: int = 300,
) -> None:
    """
    Dispersion + DOS plot (right panel).
    Saves PNG + PDF.
    """
    y_disp, ylabel = convert_units_from_cm1(freqs_cm1, unit)
    y_dos, _ = convert_units_from_cm1(dos_freq_cm1, unit)

    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    ax_dos = fig.add_subplot(gs[0, 1], sharey=ax)

    for j in range(y_disp.shape[1]):
        ax.plot(q_path, y_disp[:, j], "-", lw=2.0, color="tab:blue")
    ax.axhline(0.0, lw=0.8, color="black")

    if (emin is not None) or (emax is not None):
        lo = emin if emin is not None else float(np.nanmin(y_disp))
        hi = emax if emax is not None else float(np.nanmax(y_disp))
        ax.set_ylim(lo, hi)
    ax.set_xlim(float(np.min(q_path)), float(np.max(q_path)))

    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(full_title, fontsize=18)
    ax.tick_params(axis="y", labelsize=12)

    if labels and tick_idx and len(labels) == len(tick_idx):
        tick_pos = [q_path[i] for i in tick_idx]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(labels, fontsize=18)
        for x in tick_pos:
            ax.axvline(x, linestyle="--", linewidth=0.7, color="gray", alpha=0.6)
    else:
        ax.set_xticks([])

    ax.grid(alpha=0.3)

    # DOS panel
    ax_dos.plot(dos_val, y_dos, "-", lw=1.8, color="tab:blue")
    ax_dos.fill_betweenx(y_dos, 0.0, dos_val, alpha=0.20, color="tab:blue")
    ax_dos.set_xlabel("DOS", fontsize=12)
    ax_dos.grid(alpha=0.3)
    ax_dos.tick_params(axis="x", labelsize=10)

    plt.setp(ax_dos.get_yticklabels(), visible=False)
    ax_dos.tick_params(axis="y", length=0)

    fig.tight_layout()

    out_png = f"{out_base}.png"
    out_pdf = f"{out_base}.pdf"
    fig.savefig(out_png, dpi=dpi, transparent=True)
    fig.savefig(out_pdf, transparent=True)
    plt.close(fig)

    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")


# ==========================================================
# Main
# ==========================================================
def main() -> None:
    args = parse_args()

    freq_path, dos_path, qpath_path, prefix = resolve_inputs(args.freq, args.dos, args.qpath)
    out_base = sanitize_out_basename(args.out)

    # Load dispersion
    ext = freq_path.suffix.lower()
    if ext == ".gp":
        q_path, freqs_cm1 = load_freq_gp(freq_path)
    else:
        with open(freq_path, "r", encoding="utf-8", errors="ignore") as f:
            first = f.readline().strip()
        if first.lower().startswith("&plot"):
            q_path, freqs_cm1 = load_freq_raw(freq_path)
        else:
            q_path, freqs_cm1 = load_freq_table(freq_path)

    # HS labels (optional)
    labels: List[str] = []
    tick_idx: List[int] = []
    if qpath_path is not None:
        labels, tick_idx = parse_qpath_labels(qpath_path, len(q_path))

    # Title
    full_title = f"{args.title}"

    # DOS (optional): if missing or unreadable -> dispersion only
    if dos_path is not None:
        try:
            dos_f_cm1, dos_val = load_dos_2col(dos_path)
            plot_dispersion_with_dos(
                q_path=q_path,
                freqs_cm1=freqs_cm1,
                dos_freq_cm1=dos_f_cm1,
                dos_val=dos_val,
                out_base=out_base,
                emin=args.emin,
                emax=args.emax,
                full_title=full_title,
                unit=args.unit,
                labels=labels,
                tick_idx=tick_idx,
                dpi=args.dpi,
            )
        except Exception as e:
            print(f"[WARN] Failed to load/plot DOS from {dos_path}. Falling back to dispersion-only. ({e})")
            plot_dispersion_only(
                q_path=q_path,
                freqs_cm1=freqs_cm1,
                out_base=out_base,
                emin=args.emin,
                emax=args.emax,
                full_title=full_title,
                unit=args.unit,
                labels=labels,
                tick_idx=tick_idx,
                dpi=args.dpi,
            )
    else:
        plot_dispersion_only(
            q_path=q_path,
            freqs_cm1=freqs_cm1,
            out_base=out_base,
            emin=args.emin,
            emax=args.emax,
            full_title=full_title,
            unit=args.unit,
            labels=labels,
            tick_idx=tick_idx,
            dpi=args.dpi,
        )

    # Print resolved inputs (debug-friendly)
    print(f"[INFO] Dispersion: {freq_path}")
    print(f"[INFO] DOS      : {dos_path if dos_path is not None else 'None'}")
    print(f"[INFO] Q-path   : {qpath_path if qpath_path is not None else 'None'}")
    if prefix:
        print(f"[INFO] PREFIX   : {prefix}")


if __name__ == "__main__":
    main()
