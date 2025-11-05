#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI-only phonon dispersion plotter for Quantum ESPRESSO
-------------------------------------------------------
- Headless (Matplotlib 'Agg'), saves PNG only.
- Supports *.freq.gp (gnuplot table) and raw *.freq (&plot header).
- Parses high-symmetry (HS) point labels from QE-style q-path files.
- Draws vertical dashed lines and HS tick labels (Γ, M, K, etc).
- Styling mimics qebands.py aesthetics (clean, thicker lines, labeled HS).

QE refs:
- matdyn.x outputs frequencies for given q-list (.freq/.freq.gp).
- Standard workflow: ph.x -> q2r.x -> matdyn.x.
"""

import os, re, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- unit conversions ----
CM1_PER_THZ = 33.356
CM1_PER_MEV = 8.066

# ==========================================================
# Argument parsing
# ==========================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="Headless phonon dispersion plotter for QE (.freq.gp or .freq)."
    )
    ap.add_argument("--freq", required=True,
                    help="Path to .freq.gp (preferred) or raw .freq file.")
    ap.add_argument("--qpath", default=None,
                    help="Optional q-path file (band form) for HS labels.")
    ap.add_argument("--out", default="phonon_dispersion.png",
                    help="Output PNG filename.")
    ap.add_argument("--emin", type=float, default=None, help="Y-axis min.")
    ap.add_argument("--emax", type=float, default=None, help="Y-axis max.")
    ap.add_argument("--dpi", type=int, default=300, help="PNG DPI.")
    ap.add_argument("--title", default="Phonon dispersion", help="Figure title.")
    ap.add_argument("--unit", choices=["cm-1", "thz", "mev"], default="cm-1",
                    help="Output unit for frequencies.")
    return ap.parse_args()

# ==========================================================
# Helpers
# ==========================================================
def _is_numeric_line(s: str) -> bool:
    s = s.strip()
    if not s or s.startswith("#"):
        return False
    return bool(re.match(r'^[\s]*[+\-]?\d', s))

def load_freq_gp(path_gp: str):
    data = np.loadtxt(path_gp)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad .freq.gp table: {path_gp}")
    return data[:, 0], data[:, 1:]

def load_freq_raw(path_freq: str):
    """Parse matdyn.x raw .freq file (&plot header)."""
    with open(path_freq, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]
    m = re.match(r'^\s*&plot\s+nbnd=\s*([0-9]+)\s*,\s*nks=\s*([0-9]+)\s*/', lines[0])
    if not m:
        raise ValueError(f"Not a raw .freq header: {lines[0]}")
    nbnd, nks = int(m.group(1)), int(m.group(2))
    q_list, bands = [], []
    i = 1
    while i < len(lines) and len(q_list) < nks:
        toks = lines[i].split()
        if len(toks) != 3:
            raise ValueError(f"Unexpected q-vector line at line {i+1}")
        q_list.append(toks); i += 1
        cur = []
        while i < len(lines) and len(cur) < nbnd and _is_numeric_line(lines[i]):
            cur.extend(np.fromstring(lines[i], sep=' ').tolist()); i += 1
        bands.append(cur[:nbnd])
    q_idx = np.arange(len(q_list), dtype=float)
    return q_idx, np.array(bands, dtype=float)

def parse_qpath_labels(qpath_path: str, n_points_data: int):
    """Read HS labels and infer segment tick indices."""
    if not qpath_path or not os.path.exists(qpath_path):
        return [], []
    with open(qpath_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    try:
        n_hsp = int(lines[0].split()[0])
    except Exception:
        return [], []
    labels = []
    for i in range(1, 1 + n_hsp):
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

def convert_units(freqs_cm1: np.ndarray, unit: str):
    if unit == "cm-1":
        return freqs_cm1, "Phonon frequency (cm$^{-1}$)"
    elif unit == "thz":
        return freqs_cm1 / CM1_PER_THZ, "Phonon frequency (THz)"
    elif unit == "mev":
        return freqs_cm1 / CM1_PER_MEV, "Phonon energy (meV)"
    else:
        raise ValueError(f"Unknown unit: {unit}")

# ==========================================================
# Plot
# ==========================================================
def plot_dispersion(
    q_path, freqs_cm1, out_png,
    emin, emax, title, unit, labels, tick_idx, dpi=300
):
    freqs, ylabel = convert_units(freqs_cm1, unit)
    fig, ax = plt.subplots(figsize=(8, 6))
    # Line style
    for j in range(freqs.shape[1]):
        ax.plot(q_path, freqs[:, j], "-", lw=2.0, color="tab:blue")

    # Zero line
    ax.axhline(0.0, lw=0.8, color="black")

    # Axis limits
    if (emin is not None) or (emax is not None):
        lo = emin if emin is not None else float(np.nanmin(freqs))
        hi = emax if emax is not None else float(np.nanmax(freqs))
        ax.set_ylim(lo, hi)
    ax.set_xlim(float(np.min(q_path)), float(np.max(q_path)))

    # Labels
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(f"{title}", fontsize=18)
    ax.tick_params(axis="y", labelsize=12)

    # HS tick labeling (same style as qebands)
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
    fig.savefig(out_png, dpi=dpi, transparent=True)
    plt.close(fig)
    print(f"[OK] Saved: {out_png}")

# ==========================================================
# Main
# ==========================================================
def main():
    args = parse_args()
    ext = os.path.splitext(args.freq)[1].lower()
    if ext == ".gp":
        q_path, freqs_cm1 = load_freq_gp(args.freq)
    else:
        with open(args.freq, "r", encoding="utf-8", errors="ignore") as f:
            first = f.readline().strip()
        if first.lower().startswith("&plot"):
            q_path, freqs_cm1 = load_freq_raw(args.freq)
        else:
            data = np.loadtxt(args.freq)
            q_path, freqs_cm1 = data[:, 0], data[:, 1:]

    labels, tick_idx = ([], [])
    if args.qpath:
        labels, tick_idx = parse_qpath_labels(args.qpath, len(q_path))

    # Append unit to title
    unit_tag = {"cm-1": "(cm$^{-1}$)", "thz": "(THz)", "mev": "(meV)"}[args.unit]
    full_title = f"{args.title} {unit_tag}"

    plot_dispersion(q_path, freqs_cm1, args.out,
                    args.emin, args.emax, full_title,
                    args.unit, labels, tick_idx, dpi=args.dpi)

if __name__ == "__main__":
    main()

