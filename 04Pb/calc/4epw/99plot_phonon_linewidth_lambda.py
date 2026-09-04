#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI-only phonon dispersion + phonon DOS plotter for Quantum ESPRESSO,
with optional EPW phonon-linewidth and mode-resolved lambda overlays.
---------------------------------------------------------------------
- Headless (Matplotlib 'Agg'), saves PNG + PDF by default.
- Dispersion supports *.freq.gp (gnuplot table) and raw *.freq (&plot header).
- Auto high-symmetry (HS) labeling from qpath.in.
- Optional phonon DOS panel on the right (width ratio 5:1).
- Optional EPW overlays: linewidth.phself.*K and lambda.phself.*K are auto-detected.
- Lambda is shown as a red shaded band around each selected phonon branch.

Search behaviour
----------------
Every input (.freq.gp / .freq, .phdos, qpath.in, linewidth.phself.*K) is looked
up in the current directory first, then in the parent directory '..'.
Missing DOS / qpath / EPW files are skipped without crashing.

Typical QE workflow:
  ph.x -> q2r.x -> matdyn.x (dispersion + DOS)
  epw.x with phonselfen = .true. and filqf = <same band path>

Examples
--------
  python plot_phonon_dispersion_lw.py
  python plot_phonon_dispersion_lw.py --freq ../pb.freq.gp
  python plot_phonon_dispersion_lw.py --bubble-max 400 --gamma-max 50
  python plot_phonon_dispersion_lw.py --no-epw --emin 0
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---- unit conversions ----
CM1_PER_THZ = 33.356
CM1_PER_MEV = 8.066

# ---- directories searched for every input file ----
SEARCH_DIRS = [Path("."), Path("..")]

_RE_LW_T = re.compile(r"linewidth\.phself\.(\d+(?:\.\d+)?)K$")
_RE_LAMBDA_T = re.compile(r"lambda\.phself\.(\d+(?:\.\d+)?)K$")


# ==========================================================
# Argument parsing
# ==========================================================
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Headless phonon dispersion + DOS plotter with EPW linewidth bubbles."
    )
    ap.add_argument(
        "--freq", default=None,
        help=("Path to .freq.gp (preferred) or raw .freq file, or a PREFIX. "
              "If omitted, the script scans '.' and '..' for a dispersion file.")
    )
    ap.add_argument(
        "--dos", default=None,
        help="Optional phonon DOS file (2 columns: freq(cm^-1) DOS). Auto-detected if omitted."
    )
    ap.add_argument(
        "--qpath", default=None,
        help="Optional q-path file for HS labels. Defaults to qpath.in in '.' or '..'."
    )
    ap.add_argument(
        "--epw-qpath-out", default="qpath_epw.kpt",
        help="EPW filqf path generated automatically from qpath.in (default: qpath_epw.kpt)."
    )
    ap.add_argument(
        "--no-write-epw-qpath", action="store_true",
        help="Do not generate the EPW filqf q-path file from qpath.in."
    )
    ap.add_argument(
        "--epw-dir", default=None,
        help="Directory holding linewidth.phself.*K. Defaults to '.' then '..'."
    )
    ap.add_argument("--no-epw", action="store_true", help="Disable all EPW linewidth/lambda overlays.")
    ap.add_argument("--no-lambda", action="store_true", help="Disable the lambda overlay only.")
    ap.add_argument("--no-dos", action="store_true", help="Disable the DOS panel.")

    ap.add_argument(
        "--out", default="linewidth_lambda",
        help="Output prefix for generated figure files."
    )
    ap.add_argument("--emin", type=float, default=None, help="Y-axis min.")
    ap.add_argument("--emax", type=float, default=None, help="Y-axis max.")
    ap.add_argument("--dpi", type=int, default=300, help="PNG DPI.")
    ap.add_argument("--title", default="", help="Figure title.")
    ap.add_argument(
        "--unit", choices=["mev", "thz", "cm-1"], default="mev",
        help="Output unit for frequencies (default: meV)."
    )

    ap.add_argument("--bands", default="all",
                    help="Which dispersion branches carry the overlay: 'all' or a "
                         "1-based list such as 3,4,5. Branch n is the n-th column "
                         "of the .freq.gp file, i.e. the n-th lowest mode.")

    # ---- overlay style ----
    ap.add_argument("--style", choices=["band", "bubble"], default="band",
                    help="EPW overlay style: shaded band omega +/- A*gamma (default) or bubbles.")

    # ---- band controls ----
    ap.add_argument("--band-split", type=float, default=None,
                    help="Energy (plot unit) splitting low/high branch groups. "
                         "Omit for a single group.")
    ap.add_argument("--band-scale-lo", type=float, default=6.0,
                    help="Magnification A for the low-energy group.")
    ap.add_argument("--band-scale-hi", type=float, default=2.0,
                    help="Magnification A for the high-energy group.")
    ap.add_argument("--band-color-lo", default="dodgerblue", help="Low-group colour.")
    ap.add_argument("--band-color-hi", default="darkorange", help="High-group colour.")
    ap.add_argument("--band-alpha", type=float, default=0.35, help="Band opacity.")
    ap.add_argument("--band-lw", type=float, default=1.0, help="Dispersion line width.")

    # ---- lambda overlay controls ----
    ap.add_argument("--lambda-scale", type=float, default=2.0,
                    help="Half-width magnification for dimensionless lambda in the current plot unit.")
    ap.add_argument("--lambda-color", default="firebrick", help="Lambda overlay colour (default: firebrick).")
    ap.add_argument("--lambda-alpha", type=float, default=0.25, help="Lambda band opacity.")

    # ---- bubble controls ----
    ap.add_argument("--bubble-max", type=float, default=350.0,
                    help="Marker area (pt^2) assigned to the reference gamma.")
    ap.add_argument("--bubble-min", type=float, default=3.0,
                    help="Minimum marker area so tiny gamma stays visible.")
    ap.add_argument("--gamma-ref", type=float, default=None,
                    help="Gamma (meV) mapped to --bubble-max. Default: 95th percentile.")
    ap.add_argument("--gamma-max", type=float, default=None,
                    help="Drop points with gamma above this value (meV).")
    ap.add_argument("--bubble-alpha", type=float, default=0.35, help="Bubble opacity.")
    ap.add_argument("--bubble-color", default="crimson", help="Bubble face colour.")
    ap.add_argument("--bubble-on", choices=["qe", "epw"], default="qe",
                    help="Place bubbles on the QE branch (default) or at the raw EPW frequency.")
    return ap.parse_args()


# ==========================================================
# File auto-detection
# ==========================================================
def _find(names: List[str]) -> Optional[Path]:
    """Return the first existing file among names, searched in SEARCH_DIRS."""
    for d in SEARCH_DIRS:
        for n in names:
            p = d / n
            if p.exists() and p.is_file():
                return p
    return None


def _glob_first(pattern: str) -> Optional[Path]:
    """Return the first file matching a glob pattern, searched in SEARCH_DIRS."""
    for d in SEARCH_DIRS:
        hits = sorted(d.glob(pattern))
        for h in hits:
            if h.is_file():
                return h
    return None


def resolve_dispersion(freq_arg: Optional[str]) -> Tuple[Path, str]:
    """Locate the dispersion file. Fatal if nothing is found."""
    if freq_arg:
        p = Path(freq_arg)
        if p.exists() and p.is_file():
            return p, ""
        prefix = freq_arg
        names = [
            f"{prefix}_phband.freq.gp", f"{prefix}_phband.freq",
            f"{prefix}.freq.gp", f"{prefix}.freq",
            f"{prefix}_dispersion.freq.gp", f"{prefix}_dispersion.freq",
            f"{prefix}_phband.freq.gp.dat", f"{prefix}_phband.freq.dat",
        ]
        found = _find(names)
        if found is None:
            raise FileNotFoundError(
                f"Could not find a dispersion file for prefix '{prefix}' in "
                f"{[str(d) for d in SEARCH_DIRS]}"
            )
        return found, prefix

    # no --freq given: scan for anything that looks like a dispersion file
    for pattern in ("*.freq.gp", "*.freq"):
        hit = _glob_first(pattern)
        if hit is not None:
            return hit, hit.name.split(".")[0]
    raise FileNotFoundError(
        f"No *.freq.gp or *.freq found in {[str(d) for d in SEARCH_DIRS]}. Use --freq."
    )


def resolve_dos(dos_arg: Optional[str], prefix: str, disabled: bool) -> Optional[Path]:
    """Locate the phonon DOS file. Non-fatal."""
    if disabled:
        return None
    if dos_arg:
        p = Path(dos_arg)
        if p.exists() and p.is_file():
            return p
        print(f"[WARN] --dos not found, skipping DOS: {dos_arg}")
        return None
    names: List[str] = []
    if prefix:
        names += [f"{prefix}_phdos", f"{prefix}.phdos",
                  f"{prefix}_phdos.dat", f"{prefix}.phdos.dat"]
    found = _find(names) if names else None
    if found is None:
        found = _glob_first("*.phdos") or _glob_first("*_phdos")
    return found


def resolve_qpath(qpath_arg: Optional[str]) -> Optional[Path]:
    """Locate the q-path file used for HS labels. Non-fatal."""
    if qpath_arg:
        p = Path(qpath_arg)
        if p.exists() and p.is_file():
            return p
        print(f"[WARN] --qpath not found, skipping HS labels: {qpath_arg}")
        return None
    return _find(["qpath.in", "qpath.dat", "band.in"])


def resolve_epw(epw_dir: Optional[str], disabled: bool) -> List[Path]:
    """Locate every linewidth.phself.*K file. Non-fatal."""
    if disabled:
        return []
    dirs = [Path(epw_dir)] if epw_dir else SEARCH_DIRS
    for d in dirs:
        if not d.is_dir():
            continue
        hits = sorted(p for p in d.glob("linewidth.phself.*K")
                      if p.is_file() and _RE_LW_T.search(p.name))
        if hits:
            return hits
    return []


def resolve_lambda(epw_dir: Optional[str], disabled: bool) -> List[Path]:
    """Locate every lambda.phself.*K file. Non-fatal."""
    if disabled:
        return []
    dirs = [Path(epw_dir)] if epw_dir else SEARCH_DIRS
    for d in dirs:
        if not d.is_dir():
            continue
        hits = sorted(p for p in d.glob("lambda.phself.*K")
                      if p.is_file() and _RE_LAMBDA_T.search(p.name))
        if hits:
            return hits
    return []


def _temperature_from_path(path: Path, kind: str) -> Optional[float]:
    rx = _RE_LW_T if kind == "linewidth" else _RE_LAMBDA_T
    m = rx.search(path.name)
    return float(m.group(1)) if m else None


def sanitize_out_basename(out_arg: str) -> str:
    """Strip a .png/.pdf extension if the user supplied one."""
    p = Path(out_arg)
    if p.suffix.lower() in [".png", ".pdf"]:
        return str(p.with_suffix(""))
    return str(p)


# ==========================================================
# Parsers
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
            raise ValueError(f"Unexpected q-vector line at line {i + 1}")
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
    """Load a generic whitespace table with columns q f1 f2 ..."""
    data = np.loadtxt(str(path_any))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad frequency table: {path_any}")
    return data[:, 0], data[:, 1:]


def load_dos_2col(path_dos: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load a DOS file as two columns: freq(cm^-1) dos."""
    data = np.loadtxt(str(path_dos))
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Bad DOS table (need >=2 cols): {path_dos}")
    f_cm1 = data[:, 0].astype(float)
    dos = data[:, 1].astype(float)
    idx = np.argsort(f_cm1)
    return f_cm1[idx], dos[idx]


def load_linewidth_phself(path: Path) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read linewidth.phself.*K.

    Columns: q_index  mode_index  omega(meV)  gamma(meV)

    Returns T_K, q indices (nq,), omega (nq, nmodes) in meV, gamma (nq, nmodes) in meV.
    """
    m = _RE_LW_T.search(path.name)
    T_K = float(m.group(1)) if m else 0.0

    rows: List[Tuple[int, int, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if len(tok) < 4:
                continue
            try:
                rows.append((int(float(tok[0])), int(float(tok[1])),
                             float(tok[2]), float(tok[3])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric rows in {path.name}")

    qa = np.array([r[0] for r in rows], dtype=int)
    ma = np.array([r[1] for r in rows], dtype=int)
    wa = np.array([r[2] for r in rows], dtype=float)
    ga = np.array([r[3] for r in rows], dtype=float)

    q_uni, m_uni = np.unique(qa), np.unique(ma)
    qp = {v: i for i, v in enumerate(q_uni)}
    mp = {v: i for i, v in enumerate(m_uni)}

    omega = np.full((len(q_uni), len(m_uni)), np.nan)
    gamma = np.full((len(q_uni), len(m_uni)), np.nan)
    for qi, mi, w, g in zip(qa, ma, wa, ga):
        omega[qp[qi], mp[mi]] = w
        gamma[qp[qi], mp[mi]] = g
    return T_K, q_uni, omega, gamma


def load_lambda_phself(path: Path) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Read lambda.phself.*K. EPW format: q_index, lambda(mode 1), ..., lambda(mode N).

    Returns T_K, q indices (nq,), lambda (nq, nmodes), dimensionless.
    """
    m = _RE_LAMBDA_T.search(path.name)
    T_K = float(m.group(1)) if m else 0.0
    rows: List[List[float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            try:
                vals = [float(x) for x in s.split()]
            except ValueError:
                continue
            if len(vals) >= 2:
                rows.append(vals)
    if not rows:
        raise ValueError(f"No numeric rows in {path.name}")
    ncol = min(len(r) for r in rows)
    arr = np.array([r[:ncol] for r in rows], dtype=float)
    q_idx = arr[:, 0].astype(int)
    lam = arr[:, 1:]
    return T_K, q_idx, lam


def _clean_qpath_lines(qpath_path: Path) -> List[str]:
    """Remove blanks, comments, and Python-style triple-quoted comment blocks."""
    out: List[str] = []
    in_triple = False
    triple_token = None
    double_triple = chr(34) * 3
    single_triple = chr(39) * 3
    with qpath_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue

            if s.startswith(double_triple) or s.startswith(single_triple):
                tok = s[:3]
                if not in_triple:
                    in_triple = True
                    triple_token = tok
                    if len(s) > 3 and s.endswith(tok):
                        in_triple = False
                        triple_token = None
                    continue
                if tok == triple_token:
                    in_triple = False
                    triple_token = None
                    continue

            if in_triple:
                continue

            s = re.split(r"[#!]", s, maxsplit=1)[0].strip()
            if s:
                out.append(s)
    return out


def parse_qpath_full(qpath_path: Path) -> Tuple[List[np.ndarray], List[int], List[str]]:
    """Parse dynamic qpath.in rows: qx qy qz N_segment label."""
    lines = _clean_qpath_lines(qpath_path)
    start = None
    n_hsp = None
    for i, line in enumerate(lines):
        toks = line.split()
        if len(toks) == 1 and re.fullmatch(r"[+-]?\d+", toks[0]):
            start = i
            n_hsp = int(toks[0])
            break
    if start is None or n_hsp is None:
        raise ValueError("Could not find number of high-symmetry points")
    if n_hsp < 2:
        raise ValueError("At least two high-symmetry points are required")

    hs_lines = lines[start + 1:start + 1 + n_hsp]
    if len(hs_lines) != n_hsp:
        raise ValueError(f"Expected {n_hsp} HS rows but found {len(hs_lines)}")

    points: List[np.ndarray] = []
    counts: List[int] = []
    labels: List[str] = []
    for irow, line in enumerate(hs_lines, start=1):
        toks = line.split()
        if len(toks) < 4:
            raise ValueError(f"Bad HS row {irow}: {line}")
        q = np.array([float(toks[0]), float(toks[1]), float(toks[2])], dtype=float)
        nseg = int(float(toks[3]))
        if nseg <= 0:
            raise ValueError(f"Segment count must be positive in HS row {irow}: {line}")
        label = toks[4] if len(toks) >= 5 else f"P{irow}"
        if re.fullmatch(r"(?i)(g|gamma)", label):
            label = r"$\Gamma$"
        points.append(q)
        counts.append(nseg)
        labels.append(label)
    return points, counts, labels


def build_epw_qpath(points: List[np.ndarray], counts: List[int]) -> np.ndarray:
    """Build q points with shared HS endpoints included only once."""
    q_all: List[np.ndarray] = []
    for i in range(len(points) - 1):
        seg = np.linspace(points[i], points[i + 1], counts[i] + 1)
        if i > 0:
            seg = seg[1:]
        q_all.extend(seg)
    return np.asarray(q_all, dtype=float)


def qpath_hs_indices(counts: List[int], n_hs: int) -> List[int]:
    idx = [0]
    running = 0
    for i in range(n_hs - 1):
        running += counts[i]
        idx.append(running)
    return idx


def write_epw_qpath(qpath_path: Path, output_path: Path, weight: float = 1.0) -> Tuple[int, List[str], List[int]]:
    """Generate an EPW filqf file from qpath.in without interrupting plotting."""
    points, counts, labels = parse_qpath_full(qpath_path)
    qpts = build_epw_qpath(points, counts)
    indices = qpath_hs_indices(counts, len(points))
    with output_path.open("w", encoding="utf-8") as f:
        f.write(f"{len(qpts)} crystal\n")
        for q in qpts:
            f.write(f"{q[0]: .10f}  {q[1]: .10f}  {q[2]: .10f}  {weight:.10f}\n")
    print(f"[OK] Wrote EPW q-path: {output_path}")
    print(f"[INFO] EPW q points: {len(qpts)}")
    print("[INFO] HS points:")
    for label, idx, q in zip(labels, indices, points):
        plain_label = "G" if label == r"$\Gamma$" else label
        print(f"  {plain_label:>4s}  index={idx:5d}  q=({q[0]:.6f}, {q[1]:.6f}, {q[2]:.6f})")
    return len(qpts), labels, indices


def parse_qpath_labels(qpath_path: Path, n_points_data: int) -> Tuple[List[str], List[int]]:
    """Read HS labels and use qpath.in segment counts for exact tick positions."""
    try:
        points, counts, labels = parse_qpath_full(qpath_path)
        idx = qpath_hs_indices(counts, len(points))
        n_qpath = idx[-1] + 1
        if n_qpath != n_points_data:
            print(f"[WARN] qpath.in implies {n_qpath} dispersion points but loaded dispersion has "
                  f"{n_points_data}. Rescaling HS tick positions to the loaded path.")
            if n_qpath > 1:
                scale = (n_points_data - 1) / (n_qpath - 1)
                idx = [int(round(i * scale)) for i in idx]
        return labels, idx
    except Exception as e:
        print(f"[WARN] Failed to parse qpath labels from {qpath_path}. Skipping. ({e})")
        return [], []


def convert_units_from_cm1(y_cm1: np.ndarray, unit: str) -> Tuple[np.ndarray, str]:
    """Convert a frequency axis given in cm^-1 into the requested unit."""
    if unit == "cm-1":
        return y_cm1, "Frequency (cm$^{-1}$)"
    if unit == "thz":
        return y_cm1 / CM1_PER_THZ, "Frequency (THz)"
    if unit == "mev":
        return y_cm1 / CM1_PER_MEV, "Energy (meV)"
    raise ValueError(f"Unknown unit: {unit}")


def mev_to_unit(y_mev: np.ndarray, unit: str) -> np.ndarray:
    """Convert a meV array into the requested plotting unit."""
    out, _ = convert_units_from_cm1(np.asarray(y_mev, dtype=float) * CM1_PER_MEV, unit)
    return out


# ==========================================================
# EPW overlay
# ==========================================================
def map_epw_to_qpath(q_idx: np.ndarray, q_path: np.ndarray) -> np.ndarray:
    """
    Map EPW q indices onto the QE dispersion abscissa.

    When the two grids have the same length the mapping is one to one.
    Otherwise the indices are stretched linearly onto the path and a warning
    is printed, because that mapping is only meaningful if both calculations
    followed the same q path.
    """
    nq_epw, nq_disp = len(q_idx), len(q_path)
    if nq_epw == nq_disp:
        return q_path.astype(float)
    print(f"[WARN] EPW has {nq_epw} q points but the dispersion has {nq_disp}. "
          f"Falling back to a linear stretch. Check that filqf used the same path.")
    lo, hi = float(np.min(q_path)), float(np.max(q_path))
    span = float(q_idx.max() - q_idx.min())
    if span <= 0:
        return np.full(nq_epw, lo)
    return lo + (q_idx.astype(float) - q_idx.min()) / span * (hi - lo)


def build_bubbles(q_idx: np.ndarray, omega_mev: np.ndarray, gamma_mev: np.ndarray,
                  q_path: np.ndarray, unit: str, gamma_max: Optional[float]
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Flatten the EPW table into scatter-ready arrays.

    Returns x, y (plot unit), gamma (meV), and the number of dropped points.
    """
    x_map = map_epw_to_qpath(q_idx, q_path)
    nq, nm = omega_mev.shape

    x = np.repeat(x_map[:, None], nm, axis=1).ravel()
    w = omega_mev.ravel()
    g = np.abs(gamma_mev.ravel())

    good = np.isfinite(w) & np.isfinite(g) & (g > 0.0)
    n_drop = 0
    if gamma_max is not None:
        over = good & (g > gamma_max)
        n_drop = int(over.sum())
        good = good & ~over

    return x[good], mev_to_unit(w[good], unit), g[good], n_drop


def parse_band_selection(spec: str, n_avail: int) -> List[int]:
    """Turn a --bands string into a list of 0-based branch indices."""
    if not spec or spec.strip().lower() == "all":
        return list(range(n_avail))
    out: List[int] = []
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok and not tok.startswith("-"):
            a, b = tok.split("-", 1)
            out.extend(range(int(a) - 1, int(b)))
        else:
            out.append(int(tok) - 1)
    sel = sorted({j for j in out if 0 <= j < n_avail})
    if not sel:
        raise ValueError(f"--bands '{spec}' selected nothing (available: 1-{n_avail})")
    return sel


def gamma_on_qpath(q_idx: np.ndarray, gamma_mev: np.ndarray,
                   q_path: np.ndarray, unit: str) -> np.ndarray:
    """
    Put gamma on the dispersion abscissa, converted to the plotting unit.

    Returns an array shaped (len(q_path), nmodes). Widths are linear in the
    unit conversion, so no offset handling is needed.
    """
    g = np.abs(np.nan_to_num(gamma_mev, nan=0.0))
    g = mev_to_unit(g, unit)
    nq_epw, nm = g.shape
    if nq_epw == len(q_path):
        return g
    print(f"[WARN] EPW has {nq_epw} q points but the dispersion has {len(q_path)}. "
          f"Interpolating gamma linearly. Check that filqf used the same path.")
    x_src = np.linspace(float(q_path.min()), float(q_path.max()), nq_epw)
    out = np.empty((len(q_path), nm))
    for j in range(nm):
        out[:, j] = np.interp(q_path, x_src, g[:, j])
    return out


def draw_bands(ax, q_path: np.ndarray, y_disp: np.ndarray,
               g_disp: Optional[np.ndarray], sel: List[int], args) -> List[Patch]:
    """
    Draw every branch as a line and shade omega +/- A * gamma on the selected ones.

    Selected branches are split into a low and a high group by their mean energy
    so the two regions can use different magnifications.
    """
    nb_disp = y_disp.shape[1]

    for j in range(nb_disp):
        if j not in sel:
            ax.plot(q_path, y_disp[:, j], "-", lw=args.band_lw,
                    color="0.65", zorder=3)

    if g_disp is None:
        for j in sel:
            ax.plot(q_path, y_disp[:, j], "-", lw=args.band_lw,
                    color=args.band_color_lo, zorder=4)
        return []

    if args.band_split is None:
        groups = [(sel, args.band_scale_lo, args.band_color_lo)]
    else:
        lo = [j for j in sel if np.nanmean(y_disp[:, j]) < args.band_split]
        hi = [j for j in sel if j not in lo]
        groups = [(lo, args.band_scale_lo, args.band_color_lo),
                  (hi, args.band_scale_hi, args.band_color_hi)]

    handles: List[Patch] = []
    for idx, A, col in groups:
        if not idx:
            continue
        for j in idx:
            w = y_disp[:, j]
            half = A * g_disp[:, j]
            ax.fill_between(q_path, w - half, w + half, color=col,
                            alpha=args.band_alpha, linewidth=0, zorder=4)
            ax.plot(q_path, w, "-", lw=args.band_lw, color=col, zorder=5)
        handles.append(Patch(facecolor=col, alpha=args.band_alpha, edgecolor=col,
                             label=rf"$\omega \pm {A:g}\,\gamma$"))
    return handles


def lambda_on_qpath(q_idx: np.ndarray, lambda_vals: np.ndarray,
                    q_path: np.ndarray) -> np.ndarray:
    """Map dimensionless mode-resolved lambda onto the dispersion q grid."""
    lam = np.abs(np.nan_to_num(lambda_vals, nan=0.0))
    nq_epw, nm = lam.shape
    if nq_epw == len(q_path):
        return lam
    print(f"[WARN] lambda has {nq_epw} q points but the dispersion has {len(q_path)}. "
          f"Interpolating lambda linearly. Check that filqf used the same path.")
    x_src = np.linspace(float(q_path.min()), float(q_path.max()), nq_epw)
    out = np.empty((len(q_path), nm))
    for j in range(nm):
        out[:, j] = np.interp(q_path, x_src, lam[:, j])
    return out


def draw_lambda_bands(ax, q_path: np.ndarray, y_disp: np.ndarray,
                      lambda_disp: Optional[np.ndarray], sel: List[int], args) -> List[Patch]:
    """Draw lambda as red shaded omega +/- scale*lambda bands on selected branches."""
    if lambda_disp is None:
        return []
    handles: List[Patch] = []
    n_avail = min(y_disp.shape[1], lambda_disp.shape[1])
    used = [j for j in sel if j < n_avail]
    for j in used:
        w = y_disp[:, j]
        half = args.lambda_scale * lambda_disp[:, j]
        ax.fill_between(q_path, w - half, w + half,
                        color=args.lambda_color, alpha=args.lambda_alpha,
                        linewidth=0, zorder=5.5)
        ax.plot(q_path, w, "-", lw=args.band_lw,
                color=args.lambda_color, alpha=0.85, zorder=6)
    if used:
        handles.append(Patch(facecolor=args.lambda_color, alpha=args.lambda_alpha,
                             edgecolor=args.lambda_color,
                             label=rf"$\omega \pm {args.lambda_scale:g}\,\lambda$"))
    return handles


def draw_bubbles(ax, x: np.ndarray, y: np.ndarray, g: np.ndarray, args,
                 unit: str) -> float:
    """Draw the linewidth bubbles and return the reference gamma used."""
    if x.size == 0:
        return 0.0
    g_ref = args.gamma_ref if args.gamma_ref else float(np.percentile(g, 95))
    g_ref = max(g_ref, 1e-12)

    sizes = args.bubble_min + (args.bubble_max - args.bubble_min) * (g / g_ref)
    sizes = np.clip(sizes, args.bubble_min, args.bubble_max * 1.6)

    ax.scatter(x, y, s=sizes, facecolor=args.bubble_color, edgecolor="none",
               alpha=args.bubble_alpha, zorder=4)
    return g_ref


def bubble_legend(ax, g_ref: float, args) -> None:
    """Add a size reference legend for the bubbles."""
    if g_ref <= 0:
        return
    fracs = [0.25, 0.5, 1.0]
    handles = []
    for f in fracs:
        area = args.bubble_min + (args.bubble_max - args.bubble_min) * f
        handles.append(
            Line2D([], [], linestyle="none", marker="o",
                   markersize=np.sqrt(area), alpha=args.bubble_alpha,
                   markerfacecolor=args.bubble_color, markeredgecolor="none",
                   label=f"{f * g_ref:.3g}")
        )
    leg = ax.legend(handles=handles, title="gamma (meV)", loc="upper right",
                    fontsize=9, title_fontsize=9, labelspacing=1.4,
                    borderpad=0.9, handletextpad=1.6, framealpha=0.85)
    leg.set_zorder(10)


# ==========================================================
# Plot
# ==========================================================
def apply_common_axis_style(ax, q_path: np.ndarray, y_disp: np.ndarray,
                            labels: List[str], tick_idx: List[int],
                            ylabel: str, title: str, args,
                            show_xticklabels: bool = True,
                            show_ylabel: bool = True,
                            show_yticklabels: bool = True) -> None:
    """Apply the common cosmetics for each dispersion-style subplot."""
    ax.axhline(0.0, lw=0.8, color="black", zorder=2)

    if (args.emin is not None) or (args.emax is not None):
        lo = args.emin if args.emin is not None else float(np.nanmin(y_disp))
        hi = args.emax if args.emax is not None else float(np.nanmax(y_disp))
        ax.set_ylim(lo, hi)
    ax.set_xlim(float(np.min(q_path)), float(np.max(q_path)))

    ax.set_xlabel("")
    ax.set_title(title, fontsize=13)
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=13)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=11)
    if not show_yticklabels:
        ax.tick_params(axis="y", labelleft=False)

    if labels and tick_idx and len(labels) == len(tick_idx):
        tick_pos = [q_path[i] for i in tick_idx]
        ax.set_xticks(tick_pos)
        if show_xticklabels:
            ax.set_xticklabels(labels, fontsize=14)
        else:
            ax.set_xticklabels([])
        for xv in tick_pos:
            ax.axvline(xv, linestyle="--", linewidth=0.7, color="gray", alpha=0.6, zorder=1)
    else:
        ax.set_xticks([])

    ax.grid(alpha=0.3, zorder=0)


def draw_dispersion_base(ax, q_path: np.ndarray, y_disp: np.ndarray,
                         sel: List[int], args, selected_color: str = "tab:blue") -> None:
    """Draw the plain phonon branches before adding an overlay."""
    for j in range(y_disp.shape[1]):
        col = selected_color if j in sel else "0.65"
        ax.plot(q_path, y_disp[:, j], "-", lw=args.band_lw, color=col, zorder=3)


def draw_dos_panel(ax_dos, dos: Tuple[np.ndarray, np.ndarray], unit: str) -> None:
    """Draw the single shared DOS panel."""
    dos_f_cm1, dos_val = dos
    y_dos, _ = convert_units_from_cm1(dos_f_cm1, unit)
    ax_dos.plot(dos_val, y_dos, "-", lw=1.8, color="tab:blue")
    ax_dos.fill_betweenx(y_dos, 0.0, dos_val, alpha=0.20, color="tab:blue")
    ax_dos.set_xlabel("DOS", fontsize=11)
    ax_dos.grid(alpha=0.3)
    ax_dos.tick_params(axis="x", labelsize=10)
    plt.setp(ax_dos.get_yticklabels(), visible=False)
    ax_dos.tick_params(axis="y", length=0)


def annotate_missing_overlay(ax, message: str) -> None:
    """Put a light message inside an axis when a given overlay is missing."""
    ax.text(0.5, 0.06, message, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=10, color="0.35",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.85", alpha=0.85))


def prepare_overlay_data(q_path: np.ndarray, unit: str,
                         epw: Optional[Tuple[float, np.ndarray, np.ndarray, np.ndarray]],
                         lambda_epw: Optional[Tuple[float, np.ndarray, np.ndarray]],
                         args):
    """Prepare temperature, linewidth, and lambda arrays on the dispersion q grid."""
    T_K = None
    q_idx = omega_mev = gamma_mev = None
    g_disp = None
    lambda_disp = None

    if epw is not None:
        T_K, q_idx, omega_mev, gamma_mev = epw
        if args.gamma_max is not None:
            gamma_mev = np.where(np.abs(gamma_mev) > args.gamma_max, 0.0, gamma_mev)
        g_disp = gamma_on_qpath(q_idx, gamma_mev, q_path, unit)

    if lambda_epw is not None:
        T_lam, q_idx_lam, lambda_vals = lambda_epw
        if T_K is None:
            T_K = T_lam
        lambda_disp = lambda_on_qpath(q_idx_lam, lambda_vals, q_path)

    return T_K, q_idx, omega_mev, gamma_mev, g_disp, lambda_disp


def format_temperature_tag(T: float) -> str:
    return f"T{T:08.3f}K"


def make_chunk_out_base(out_prefix: str, temp_entries: List[dict]) -> str:
    """Build a readable output basename from the temperatures contained in a page."""
    tags = [format_temperature_tag(entry["T"]) for entry in temp_entries]
    return f"{out_prefix}_" + "_".join(tags)


def make_plain_plot(q_path: np.ndarray, freqs_cm1: np.ndarray,
                    dos: Optional[Tuple[np.ndarray, np.ndarray]],
                    labels: List[str], tick_idx: List[int],
                    out_base: str, args) -> None:
    """Fallback plot when no linewidth/lambda data are present."""
    unit = args.unit
    y_disp, ylabel = convert_units_from_cm1(freqs_cm1, unit)

    if dos is not None:
        fig = plt.figure(figsize=(10, 4.4))
        gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        ax_dos = fig.add_subplot(gs[0, 1], sharey=ax)
    else:
        fig, ax = plt.subplots(figsize=(8, 4.4))
        ax_dos = None

    sel = list(range(y_disp.shape[1]))
    draw_dispersion_base(ax, q_path, y_disp, sel, args, selected_color="tab:blue")
    title = args.title if args.title else "Phonon dispersion"
    apply_common_axis_style(ax, q_path, y_disp, labels, tick_idx, ylabel, title, args, show_xticklabels=True)

    if ax_dos is not None and dos is not None:
        draw_dos_panel(ax_dos, dos, unit)

    fig.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.12)
    out_png, out_pdf = f"{out_base}.png", f"{out_base}.pdf"
    fig.savefig(out_png, dpi=args.dpi, transparent=True)
    fig.savefig(out_pdf, transparent=True)
    plt.close(fig)
    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")


def make_grouped_plot(q_path: np.ndarray, freqs_cm1: np.ndarray,
                      dos: Optional[Tuple[np.ndarray, np.ndarray]],
                      temp_entries: List[dict],
                      labels: List[str], tick_idx: List[int],
                      out_base: str, args, total_temps: int) -> None:
    """Plot up to two temperatures per file, with linewidth left, lambda middle, DOS right."""
    unit = args.unit
    y_disp, ylabel = convert_units_from_cm1(freqs_cm1, unit)

    n_used = len(temp_entries)
    nrows = 2 if total_temps > 2 else n_used
    nrows = max(1, nrows)
    has_dos = dos is not None

    fig_w = 13 if has_dos else 11
    fig_h = 4.3 * nrows
    ncols = 3 if has_dos else 2
    width_ratios = [4.7, 4.7, 1.3] if has_dos else [1, 1]
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(nrows, ncols, width_ratios=width_ratios, hspace=0.20, wspace=0.08)

    row_axes = []
    for i in range(nrows):
        ax_lw = fig.add_subplot(gs[i, 0])
        ax_lam = fig.add_subplot(gs[i, 1], sharex=ax_lw, sharey=ax_lw)
        ax_dos = fig.add_subplot(gs[i, 2], sharey=ax_lw) if has_dos else None
        row_axes.append((ax_lw, ax_lam, ax_dos))

    last_filled_row = n_used - 1

    for i in range(nrows):
        ax_lw, ax_lam, ax_dos = row_axes[i]

        if i >= n_used:
            ax_lw.axis("off")
            ax_lam.axis("off")
            if ax_dos is not None:
                ax_dos.axis("off")
            continue

        entry = temp_entries[i]
        T_nom = entry["T"]
        epw = entry["epw"]
        lambda_epw = entry["lambda"]
        T_K, q_idx, omega_mev, gamma_mev, g_disp, lambda_disp = prepare_overlay_data(
            q_path, unit, epw, lambda_epw, args
        )
        T_show = T_K if T_K is not None else T_nom

        n_avail = y_disp.shape[1]
        if g_disp is not None:
            n_avail = min(n_avail, g_disp.shape[1])
        if lambda_disp is not None:
            n_avail = min(n_avail, lambda_disp.shape[1])
        sel = parse_band_selection(args.bands, n_avail)
        print(f"[INFO] T = {T_show:.3f} K, overlay on branches (1-based): {[j + 1 for j in sel]}")

        show_xticklabels = (i == last_filled_row)

        # ---- left panel: linewidth
        if args.style == "band":
            band_handles = draw_bands(ax_lw, q_path, y_disp, g_disp, sel, args)
            if g_disp is None:
                annotate_missing_overlay(ax_lw, "No linewidth data")
            elif band_handles:
                ax_lw.legend(handles=band_handles, loc="upper right", fontsize=10,
                             framealpha=0.85).set_zorder(10)
        else:
            draw_dispersion_base(ax_lw, q_path, y_disp, sel, args, selected_color="tab:blue")
            if g_disp is None:
                annotate_missing_overlay(ax_lw, "No linewidth data")
            else:
                if args.bubble_on == "epw":
                    x, y, g_arr, n_drop = build_bubbles(
                        q_idx, omega_mev, gamma_mev, q_path, unit, args.gamma_max
                    )
                else:
                    xs, ys, gs_ = [], [], []
                    for j in sel:
                        gj = g_disp[:, j]
                        good = np.isfinite(gj) & (gj > 0.0)
                        xs.append(q_path[good])
                        ys.append(y_disp[good, j])
                        gs_.append(gj[good])
                    x = np.concatenate(xs) if xs else np.array([])
                    y = np.concatenate(ys) if ys else np.array([])
                    g_arr = np.concatenate(gs_) if gs_ else np.array([])
                g_ref = draw_bubbles(ax_lw, x, y, g_arr, args, unit)
                print(f"[INFO] linewidth bubbles drawn = {x.size}, gamma_ref = {g_ref:.6g}")
                bubble_legend(ax_lw, g_ref, args)

        lw_title = f"Linewidth, T = {T_show:.3f} K"
        if args.title:
            lw_title = f"{args.title} — linewidth, T = {T_show:.3f} K"
        apply_common_axis_style(ax_lw, q_path, y_disp, labels, tick_idx,
                                ylabel, lw_title, args,
                                show_xticklabels=show_xticklabels,
                                show_ylabel=True,
                                show_yticklabels=True)

        # ---- middle panel: lambda
        draw_dispersion_base(ax_lam, q_path, y_disp, sel, args, selected_color=args.lambda_color)
        lambda_handles = draw_lambda_bands(ax_lam, q_path, y_disp, lambda_disp, sel, args)
        if lambda_disp is None:
            annotate_missing_overlay(ax_lam, "No lambda data")
        elif lambda_handles:
            ax_lam.legend(handles=lambda_handles, loc="upper right", fontsize=10,
                          framealpha=0.85).set_zorder(10)

        lam_title = f"Lambda, T = {T_show:.3f} K"
        if args.title:
            lam_title = f"{args.title} — lambda, T = {T_show:.3f} K"
        apply_common_axis_style(ax_lam, q_path, y_disp, labels, tick_idx,
                                ylabel, lam_title, args,
                                show_xticklabels=show_xticklabels,
                                show_ylabel=False,
                                show_yticklabels=False)

        # ---- right panel: DOS, once per row and sharing the same y-scale
        if ax_dos is not None and dos is not None:
            draw_dos_panel(ax_dos, dos, unit)
            ax_dos.set_title("Phonon DOS" if i == 0 else "", fontsize=13)
            if not show_xticklabels:
                ax_dos.set_xlabel("")
                ax_dos.tick_params(axis="x", labelbottom=False)

    fig.subplots_adjust(left=0.07, right=0.98, top=0.94, bottom=0.09)
    out_png, out_pdf = f"{out_base}.png", f"{out_base}.pdf"
    fig.savefig(out_png, dpi=args.dpi, transparent=True)
    fig.savefig(out_pdf, transparent=True)
    plt.close(fig)
    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")


# ==========================================================
# Main
# ==========================================================
def main() -> None:
    args = parse_args()

    freq_path, prefix = resolve_dispersion(args.freq)
    dos_path = resolve_dos(args.dos, prefix, args.no_dos)
    qpath_path = resolve_qpath(args.qpath)
    epw_files = resolve_epw(args.epw_dir, args.no_epw)
    lambda_files = resolve_lambda(args.epw_dir, args.no_epw or args.no_lambda)
    out_prefix = sanitize_out_basename(args.out)

    # ---- generate EPW q-path, then continue the original plotting workflow
    if qpath_path is not None and not args.no_write_epw_qpath:
        try:
            write_epw_qpath(qpath_path, Path(args.epw_qpath_out))
        except Exception as e:
            print(f"[WARN] Failed to generate EPW q-path from {qpath_path}. Continuing plotting. ({e})")

    print(f"[INFO] Dispersion: {freq_path}")
    print(f"[INFO] DOS       : {dos_path if dos_path else 'None'}")
    print(f"[INFO] Q-path    : {qpath_path if qpath_path else 'None'}")
    print(f"[INFO] linewidth : {[p.name for p in epw_files] if epw_files else 'None'}")
    print(f"[INFO] lambda    : {[p.name for p in lambda_files] if lambda_files else 'None'}")

    # ---- dispersion
    if freq_path.suffix.lower() == ".gp":
        q_path, freqs_cm1 = load_freq_gp(freq_path)
    else:
        with open(freq_path, "r", encoding="utf-8", errors="ignore") as f:
            first = f.readline().strip()
        if first.lower().startswith("&plot"):
            q_path, freqs_cm1 = load_freq_raw(freq_path)
        else:
            q_path, freqs_cm1 = load_freq_table(freq_path)
    print(f"[INFO] dispersion: {len(q_path)} q points, {freqs_cm1.shape[1]} branches")

    # ---- HS labels
    labels, tick_idx = ([], [])
    if qpath_path is not None:
        labels, tick_idx = parse_qpath_labels(qpath_path, len(q_path))

    # ---- DOS
    dos = None
    if dos_path is not None:
        try:
            dos = load_dos_2col(dos_path)
        except Exception as e:
            print(f"[WARN] Failed to load DOS from {dos_path}, skipping panel. ({e})")
            dos = None

    # ---- no EPW data: fallback plain dispersion plot
    if not epw_files and not lambda_files:
        make_plain_plot(q_path, freqs_cm1, dos, labels, tick_idx,
                        f"{out_prefix}_no_epw", args)
        return

    # ---- match linewidth/lambda files by temperature
    lw_by_T = {_temperature_from_path(p, "linewidth"): p for p in epw_files}
    lam_by_T = {_temperature_from_path(p, "lambda"): p for p in lambda_files}
    temps = sorted(t for t in (set(lw_by_T) | set(lam_by_T)) if t is not None)

    temp_entries: List[dict] = []
    for T in temps:
        epw_data = None
        lambda_data = None

        p_lw = lw_by_T.get(T)
        if p_lw is not None:
            try:
                epw_data = load_linewidth_phself(p_lw)
                T_K, q_idx, omega_mev, gamma_mev = epw_data
                n_imag = int(np.sum(omega_mev < 0))
                print(f"\n[INFO] {p_lw.name}: T = {T_K:.3f} K, "
                      f"{omega_mev.shape[0]} q points, {omega_mev.shape[1]} modes, "
                      f"imaginary points = {n_imag}")
            except Exception as e:
                print(f"[WARN] Failed to read {p_lw.name}, skipping linewidth. ({e})")

        p_lam = lam_by_T.get(T)
        if p_lam is not None:
            try:
                lambda_data = load_lambda_phself(p_lam)
                T_lam, q_lam, lam = lambda_data
                print(f"[INFO] {p_lam.name}: T = {T_lam:.3f} K, "
                      f"{lam.shape[0]} q points, {lam.shape[1]} modes, "
                      f"lambda max = {np.nanmax(lam):.6g}")
            except Exception as e:
                print(f"[WARN] Failed to read {p_lam.name}, skipping lambda. ({e})")

        if epw_data is None and lambda_data is None:
            continue
        temp_entries.append({"T": T, "epw": epw_data, "lambda": lambda_data})

    if not temp_entries:
        make_plain_plot(q_path, freqs_cm1, dos, labels, tick_idx,
                        f"{out_prefix}_no_valid_epw", args)
        return

    # ---- make at most two temperature rows per file, low T first
    total_temps = len(temp_entries)
    for i in range(0, total_temps, 2):
        chunk = temp_entries[i:i + 2]
        out_base = make_chunk_out_base(out_prefix, chunk)
        make_grouped_plot(q_path, freqs_cm1, dos, chunk,
                          labels, tick_idx, out_base, args,
                          total_temps=total_temps)


if __name__ == "__main__":
    main()
