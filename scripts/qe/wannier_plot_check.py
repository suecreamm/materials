#!/usr/bin/env python3
"""
compare_index_norm_with_fermi.py

Use k-distance (from QE k-points) for DFT x-axis, so DFT and Wannier share
the same definition of "path coordinate" before normalizing to [0,1].

Default behavior:
- Normalize path to [0,1] for both DFT and Wannier.
- Resolve Ef from:
    1) --set-fermi
    2) --fermi-from FILE
    3) auto-search '*nscf*.out' then '*scf*.out'
  If none found, skip alignment silently (no crash).
- Apply SAME Ef to BOTH DFT and Wannier by default: E <- E - Ef
- Plot Wannier blocks separately (each block is typically one band).

(Added)
- Read *.labelinfo.dat to place high-symmetry labels on x-axis and draw vertical lines.

Usage:
  python compare_index_norm_with_fermi.py band.dat TiS2_monolayer_band.dat
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path
from typing import List, Tuple, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# 1) QE &plot band parser (returns k-points + energies)
# ============================================================

def parse_qe_plot_format(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse QE band file in '&plot nbnd=..., nks=... /' format.

    Returns
    -------
    kpts : np.ndarray, shape (nks, 3)
        Fractional k-points as stored in the file.
    E : np.ndarray, shape (nbnd, nks)
        Band energies (no Ef correction here).
    """
    with open(filename, "r") as f:
        header = None
        for line in f:
            if line.strip():
                header = line.strip()
                break

        if header is None or (not header.startswith("&plot")):
            raise ValueError("Not a QE '&plot' band file.")

        m_nbnd = re.search(r"nbnd\s*=\s*(\d+)", header)
        m_nks = re.search(r"nks\s*=\s*(\d+)", header)
        if m_nbnd is None or m_nks is None:
            raise ValueError("Failed to parse nbnd/nks from '&plot' header.")

        nbnd = int(m_nbnd.group(1))
        nks = int(m_nks.group(1))

        def read_floats(n: int) -> List[float]:
            vals: List[float] = []
            while len(vals) < n:
                parts = next(f).split()
                vals.extend([float(p) for p in parts])
            return vals

        kpts = np.zeros((nks, 3), dtype=float)
        E = np.zeros((nbnd, nks), dtype=float)
        for ik in range(nks):
            kpts[ik] = read_floats(3)
            E[:, ik] = read_floats(nbnd)

    return kpts, E


def kdist_normalized_from_kpts(kpts: np.ndarray) -> np.ndarray:
    """
    Build cumulative k-distance along the path and normalize to [0,1].

    Notes
    -----
    Uses Euclidean distance in the coordinates stored in the band file.
    This is usually sufficient for aligning to Wannier 'k-distance' x.
    """
    if kpts.shape[0] <= 1:
        return np.zeros(kpts.shape[0], dtype=float)

    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(dk)])
    span = float(cum.max() - cum.min())
    if span == 0.0:
        return np.zeros_like(cum)
    return (cum - float(cum.min())) / span


# ============================================================
# 1.5) labelinfo.dat parser (ADDED)
# ============================================================

def parse_labelinfo_positions(labelinfo_path: str, dft_x: np.ndarray) -> Tuple[List[float], List[str]]:
    """
    Parse QE label info file '*.labelinfo.dat' and convert to x-positions on the
    normalized DFT path.

    Expected format (per line; extra columns allowed):
      LABEL   INDEX   KDIST   KX   KY   KZ   ...

    We trust INDEX (1-based k-point index). The x-position is taken from dft_x[INDEX-1].

    Returns
    -------
    xs : list of float
        x positions in [0,1] for each label.
    labs : list of str
        corresponding labels (e.g., G, M, K, G).
    """
    xs: List[float] = []
    labs: List[str] = []

    if labelinfo_path is None:
        return xs, labs

    p = Path(labelinfo_path)
    if not p.is_file():
        return xs, labs

    nks = int(dft_x.size)
    with open(p, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue

            lab = parts[0]
            try:
                idx_1based = int(parts[1])
            except ValueError:
                continue

            if idx_1based < 1 or idx_1based > nks:
                # Out-of-range index; ignore silently (robustness)
                continue

            xs.append(float(dft_x[idx_1based - 1]))
            labs.append(lab)

    return xs, labs


# ============================================================
# 2) Wannier 2-column parser (blocks)
# ============================================================

def parse_wannier_2col_blocks(filename: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Parse Wannier90 gnuplot-style 2-column file into blocks separated by blank lines/'e'.

    Returns
    -------
    x_blocks, y_blocks : list of arrays
        Each block typically corresponds to one band.
    """
    x_blocks: List[np.ndarray] = []
    y_blocks: List[np.ndarray] = []
    cur_x: List[float] = []
    cur_y: List[float] = []

    def flush() -> None:
        nonlocal cur_x, cur_y
        if cur_x:
            x_blocks.append(np.array(cur_x, dtype=float))
            y_blocks.append(np.array(cur_y, dtype=float))
            cur_x, cur_y = [], []

    with open(filename, "r") as f:
        for line in f:
            s = line.strip()
            if (not s) or s == "e" or s.startswith("#"):
                flush()
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                cur_x.append(float(parts[0]))
                cur_y.append(float(parts[1]))
            except ValueError:
                continue
    flush()
    return x_blocks, y_blocks


def normalize_wannier_x_blocks(x_blocks: List[np.ndarray]) -> List[np.ndarray]:
    """Normalize Wannier x (k-distance) to [0,1] using global min/max across all blocks."""
    if not x_blocks:
        return x_blocks
    all_x = np.concatenate([xb for xb in x_blocks if xb.size > 0])
    xmin = float(np.nanmin(all_x))
    xmax = float(np.nanmax(all_x))
    span = xmax - xmin
    if (not np.isfinite(span)) or span == 0.0:
        return x_blocks
    return [(xb - xmin) / span for xb in x_blocks]


# ============================================================
# 3) Ef parsing & search
# ============================================================

_FERMI_RE = re.compile(
    r"the\s+Fermi\s+energy\s+is\s+([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*ev",
    re.IGNORECASE,
)

def parse_fermi_from_qe_out(path: str) -> Optional[float]:
    """Return Ef in eV if found, else None."""
    try:
        ef: Optional[float] = None
        with open(path, "r", errors="ignore") as f:
            for line in f:
                m = _FERMI_RE.search(line)
                if m:
                    ef = float(m.group(1))
        return ef
    except Exception:
        return None


def find_qe_out_for_fermi(cwd: str = ".") -> Optional[str]:
    """Search '*nscf*.out' then '*scf*.out'; return first match or None."""
    for pat in ("*nscf*.out", "*scf*.out"):
        cands = sorted(glob.glob(str(Path(cwd) / pat)))
        if cands:
            return cands[0]
    return None


def resolve_ef(no_fermi_search: bool,
               fermi_from: Optional[str],
               set_fermi: Optional[float]) -> Tuple[Optional[float], str]:
    """Resolve Ef with priority: --set-fermi > --fermi-from > auto-search (unless disabled)."""
    if set_fermi is not None:
        return float(set_fermi), "manual(--set-fermi)"

    if fermi_from is not None:
        ef = parse_fermi_from_qe_out(fermi_from)
        if ef is None:
            return None, f"fermi-from({Path(fermi_from).name}, no-match)"
        return ef, f"fermi-from({Path(fermi_from).name})"

    if no_fermi_search:
        return None, "disabled(--no-fermi-search)"

    outp = find_qe_out_for_fermi(".")
    if outp is None:
        return None, "auto(no scf/nscf out found)"
    ef = parse_fermi_from_qe_out(outp)
    if ef is None:
        return None, f"auto({Path(outp).name}, no-match)"
    return ef, f"auto({Path(outp).name})"


# ============================================================
# 4) Plot
# ============================================================

def plot_bands(dft_file: str,
               wann_file: str,
               outpng: str,
               ylim: Tuple[float, float],
               ef: Optional[float],
               ef_src: str,
               align_fermi: bool,
               wannier_ef_override: Optional[float],
               labelinfo_file: Optional[str],  # (ADDED)
               label_fontsize: int = 16,
               tick_fontsize: int = 14) -> None:
    """
    Plot DFT vs Wannier bands with normalized path, optional Ef alignment.
    """

    # ---- DFT ----
    kpts, dft_E = parse_qe_plot_format(dft_file)
    dft_x = kdist_normalized_from_kpts(kpts)

    # ---- Wannier ----
    wx_blocks, wy_blocks = parse_wannier_2col_blocks(wann_file)
    wx_blocks = normalize_wannier_x_blocks(wx_blocks)

    # ---- Apply Ef to BOTH by default ----
    if align_fermi and (ef is not None):
        dft_E = dft_E - ef
        ef_w = ef if wannier_ef_override is None else float(wannier_ef_override)
        wy_blocks = [yb - ef_w for yb in wy_blocks]

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(6, 6))

    for b in range(dft_E.shape[0]):
        label = "DFT (QE)" if b == 0 else None
        ax.plot(dft_x, dft_E[b], color="0.4", linewidth=2.0, alpha=0.6, label=label)

    for i, (xb, yb) in enumerate(zip(wx_blocks, wy_blocks)):
        label = "Wannier (MLWF)" if i == 0 else None
        ax.plot(xb, yb, "r--", linewidth=1.5, alpha=1.0, label=label)

    # ---- High-symmetry labels from *.labelinfo.dat (ADDED) ----
    xs, labs = parse_labelinfo_positions(labelinfo_file, dft_x)
    if xs:
        # vertical separators
        for x in xs:
            ax.axvline(x, color="0.75", linewidth=1.0, linestyle="-", alpha=0.9, zorder=0)

        # set x ticks to high-sym points
        ax.set_xticks(xs)
        ax.set_xticklabels(labs, fontsize=tick_fontsize)

    # Styling
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel("Normalized Path", fontsize=label_fontsize)
    ax.set_ylabel("Energy (eV)", fontsize=label_fontsize)
    ax.grid(True, linestyle=":", alpha=0.6)

    # Tick label font sizes (do both for consistency)
    # (If labelinfo is present, xticklabels already set above; this keeps y consistent.)
    ax.tick_params(axis="y", labelsize=tick_fontsize)

    ax.legend()

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    print(f"Saved to {outpng}")
    print("[info] Ef:", "None" if ef is None else f"{ef:.6f} eV", f"({ef_src})")
    print("[info] Wannier blocks plotted:", len(wx_blocks))
    if labelinfo_file:
        p = Path(labelinfo_file)
        print("[info] labelinfo:", str(p) if p.is_file() else f"{labelinfo_file} (not found)")


def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("dft_band_file")
    p.add_argument("wann_band_file")
    p.add_argument("--out", default="band_comparison.png")

    # Default ylim = [-3, 3]
    p.add_argument("--ylim", nargs=2, type=float, default=[-3, 3], metavar=("YMIN", "YMAX"))

    p.add_argument("--no-fermi-search", action="store_true")
    p.add_argument("--fermi-from", default=None)
    p.add_argument("--set-fermi", type=float, default=None)

    p.add_argument("--no-align-fermi", action="store_true",
                   help="Disable Ef alignment for both DFT and Wannier.")
    p.add_argument("--wannier-fermi", type=float, default=None,
                   help="Override Ef used for Wannier only (rare).")

    # (ADDED) labelinfo
    p.add_argument("--labelinfo", default=None,
                   help="Path to '*.labelinfo.dat' for high-symmetry x-axis labels.")

    # Font sizes
    p.add_argument("--label-fontsize", type=int, default=16,
                   help="Font size for axis labels (default: 16).")
    p.add_argument("--tick-fontsize", type=int, default=14,
                   help="Font size for tick labels (default: 14).")

    args = p.parse_args(argv)

    # (ADDED) auto labelinfo if not provided:
    # Try: <wann_band_file>.labelinfo.dat  OR  <stem>.labelinfo.dat in cwd
    labelinfo_file = args.labelinfo
    if labelinfo_file is None:
        # common naming you showed: TiS2_monolayer_band.labelinfo.dat
        stem = Path(args.wann_band_file).name
        cand1 = Path(f"{stem}.labelinfo.dat")
        cand2 = Path(f"{Path(stem).stem}.labelinfo.dat")
        if cand1.is_file():
            labelinfo_file = str(cand1)
        elif cand2.is_file():
            labelinfo_file = str(cand2)
        else:
            # also try wildcard in cwd (best-effort)
            hits = sorted(Path(".").glob("*.labelinfo.dat"))
            if hits:
                labelinfo_file = str(hits[0])

    ef, ef_src = resolve_ef(args.no_fermi_search, args.fermi_from, args.set_fermi)
    align_fermi = (not args.no_align_fermi)

    plot_bands(
        dft_file=args.dft_band_file,
        wann_file=args.wann_band_file,
        outpng=args.out,
        ylim=(args.ylim[0], args.ylim[1]),
        ef=ef,
        ef_src=ef_src,
        align_fermi=align_fermi,
        wannier_ef_override=args.wannier_fermi,
        labelinfo_file=labelinfo_file,  # (ADDED)
        label_fontsize=args.label_fontsize,
        tick_fontsize=args.tick_fontsize,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

