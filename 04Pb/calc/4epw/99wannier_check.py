#!/usr/bin/env python3
"""
wannier_plot_check.py

Band comparison (QE DFT vs Wannier90 MLWF) + PDOS/DOS panel.

Left panel (width 4):
- DFT bands from QE '&plot' band file (band.dat)
- Wannier bands from Wannier90 gnuplot-style 2-col blocks (graphene_band.dat)

Right panel (width 1):
- Total DOS (from '*pdos_tot') in red
- Projected DOS grouped by (species-instance, orbital letter), e.g. Ti#1 d, S#1 p, S#2 s
  Colors:
    - Each chemical species gets one colormap (spring, summer, winter, autumn, ...)
    - For each species-instance (e.g. S#1), orbitals (s/p/d/...) get evenly spaced colors
      sampled from that species' colormap.

Fermi alignment:
- Resolve Ef from:
    1) --set-fermi
    2) --fermi-from FILE
    3) auto-search '*nscf*.out' then '*scf*.out'
- If found and alignment enabled: subtract Ef from BOTH band energies and PDOS energy axis.

High-symmetry labels:
- Read '*_band.labelinfo.dat' (rule: <prefix>_band.labelinfo.dat)
- Use KDIST (3rd column), normalized to [0,1], to place ticks/vertical lines.

Usage:
  python wannier_plot_check.py band.dat graphene_band.dat
  python wannier_plot_check.py band.dat band.eig --epw
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# 1) QE '&plot' band parser (returns k-points + energies)
# ============================================================

def parse_qe_plot_format(filename: str) -> Tuple[np.ndarray, np.ndarray]:
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


def kdist_cumulative_raw(kpts: np.ndarray) -> np.ndarray:
    if kpts.shape[0] <= 1:
        return np.zeros(kpts.shape[0], dtype=float)
    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(dk)])
    return cum


def normalize_01(x: np.ndarray) -> np.ndarray:
    if x.size == 0:
        return x
    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    span = xmax - xmin
    if (not np.isfinite(span)) or span == 0.0:
        return np.zeros_like(x)
    return (x - xmin) / span


# ============================================================
# 1.5) labelinfo.dat parser
# ============================================================

@dataclass
class LabelInfo:
    label: str
    idx_1based: int
    kdist: float
    kfrac: Optional[np.ndarray] = None


def read_labelinfo_file(labelinfo_path: str, debug: bool = False) -> List[LabelInfo]:
    """
    Expected format: LABEL   INDEX   KDIST   KX   KY   KZ ...
    """
    p = Path(labelinfo_path)
    if not p.is_file():
        return []

    out: List[LabelInfo] = []
    with open(p, "r") as f:
        for ln, line in enumerate(f, start=1):
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue

            lab = parts[0]
            try:
                idx = int(parts[1])
                kd = float(parts[2])
            except ValueError:
                continue

            kf: Optional[np.ndarray] = None
            if len(parts) >= 6:
                try:
                    kf = np.array([float(parts[3]), float(parts[4]), float(parts[5])], dtype=float)
                except ValueError:
                    kf = None

            out.append(LabelInfo(label=lab, idx_1based=idx, kdist=kd, kfrac=kf))

    return out


def label_x_positions_from_kdist(
    labels: List[LabelInfo],
    debug: bool = False
) -> Tuple[List[float], List[str]]:
    """
    Normalize KDIST to [0,1] to match band-path normalization.
    """
    xs: List[float] = []
    labs: List[str] = []

    if not labels:
        return xs, labs

    kds = np.array([li.kdist for li in labels], dtype=float)
    kmin = float(np.nanmin(kds))
    kmax = float(np.nanmax(kds))

    if kmax == kmin:
        return [], []

    for li in labels:
        x = (li.kdist - kmin) / (kmax - kmin)
        xs.append(x)
        labs.append(li.label)
        if debug:
            print(f"[debug] label {li.label} mapped to x={x:.6f} via KDIST={li.kdist:.6f}")

    return xs, labs


# ============================================================
# 2) Wannier 2-column parser (blocks)
# ============================================================

def parse_wannier_2col_blocks(filename: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
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


def parse_epw_band_eig(filename: str) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Parse EPW band.eig (QE/plotband-compatible '&plot' format) and return
    band-by-band x/y blocks, matching the interface used for Wannier90 data.
    """
    kpts, E = parse_qe_plot_format(filename)
    x = normalize_01(kdist_cumulative_raw(kpts))
    x_blocks = [x.copy() for _ in range(E.shape[0])]
    y_blocks = [E[ib, :].copy() for ib in range(E.shape[0])]
    return x_blocks, y_blocks



def detect_path_boundaries(
    kpts: np.ndarray,
    angle_tol: float = 1.0e-3,
    min_step: float = 1.0e-10,
) -> List[int]:
    """
    Detect high-symmetry path boundaries from changes in direction of the
    k-point sequence.  This uses only collinearity in the supplied coordinate
    representation, so it does not assume an orthogonal reciprocal basis.

    'angle_tol' is a threshold on sin(angle) between consecutive steps.  It must
    stay well above the round-off level of the printed k-coordinates (band.dat
    and band.eig carry ~6 decimals), otherwise every point looks like a corner.

    Returns 0-based point indices including the first and last point.
    """
    n = kpts.shape[0]
    if n <= 1:
        return [0]

    boundaries = [0]
    prev_dir: Optional[np.ndarray] = None

    for i in range(n - 1):
        step = kpts[i + 1] - kpts[i]
        norm = float(np.linalg.norm(step))
        if norm <= min_step:
            # A repeated k-point marks an explicit break in the path.
            if boundaries[-1] != i:
                boundaries.append(i)
            prev_dir = None
            continue
        cur_dir = step / norm

        if prev_dir is not None:
            # A segment boundary is where two consecutive non-zero steps are
            # no longer collinear (parallel or antiparallel).
            cross_norm = float(np.linalg.norm(np.cross(prev_dir, cur_dir)))
            dot = float(np.dot(prev_dir, cur_dir))
            if cross_norm > angle_tol or dot < 0.0:
                if boundaries[-1] != i:
                    boundaries.append(i)
        prev_dir = cur_dir

    if boundaries[-1] != n - 1:
        boundaries.append(n - 1)
    return boundaries


def piecewise_common_x(
    kpts: np.ndarray,
    boundaries: List[int],
    anchors: np.ndarray,
) -> np.ndarray:
    """
    Map each path segment onto common anchor positions.  Distances are only
    used *inside* a straight segment; the high-symmetry endpoints themselves
    are fixed by anchors.  This is useful when QE and EPW express equivalent
    k-points in different reciprocal-coordinate conventions.
    """
    if len(boundaries) != len(anchors):
        raise ValueError("Number of path boundaries and anchor positions differ.")

    x = np.zeros(kpts.shape[0], dtype=float)

    for iseg in range(len(boundaries) - 1):
        i0 = boundaries[iseg]
        i1 = boundaries[iseg + 1]
        if i1 <= i0:
            continue

        seg = kpts[i0:i1 + 1]
        local = kdist_cumulative_raw(seg)
        if local[-1] > 0.0:
            local = local / local[-1]
        else:
            local = np.linspace(0.0, 1.0, len(seg))

        xa = float(anchors[iseg])
        xb = float(anchors[iseg + 1])
        x[i0:i1 + 1] = xa + (xb - xa) * local

    return x


def normalize_wannier_x_blocks(x_blocks: List[np.ndarray]) -> List[np.ndarray]:
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
# 2.5) k-coordinate convention handling
# ============================================================

_B_VEC_RE = re.compile(
    r"b\(([123])\)\s*=\s*\(\s*"
    r"([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)\s*\)"
)


def parse_reciprocal_axes(path: Optional[str]) -> Optional[np.ndarray]:
    """
    Read b(1), b(2), b(3) (cart. coord. in units 2 pi/alat) from a QE output.
    Returns a 3x3 array whose ROWS are b1, b2, b3, or None.
    """
    if path is None:
        return None
    try:
        rows: Dict[int, List[float]] = {}
        with open(path, "r", errors="ignore") as f:
            for line in f:
                m = _B_VEC_RE.search(line)
                if m:
                    i = int(m.group(1))
                    if i not in rows:
                        rows[i] = [
                            float(m.group(2).replace("D", "E").replace("d", "e")),
                            float(m.group(3).replace("D", "E").replace("d", "e")),
                            float(m.group(4).replace("D", "E").replace("d", "e")),
                        ]
        if len(rows) != 3:
            return None
        return np.array([rows[1], rows[2], rows[3]], dtype=float)
    except Exception:
        return None


def label_index_mismatch(kpts: np.ndarray, labels: List[LabelInfo]) -> Optional[float]:
    """
    labelinfo.dat carries both the point INDEX along the band k-path and the
    physical path length KDIST (1/Ang) of each high-symmetry point.  If a given
    k-coordinate representation is the metrically correct one, then the
    normalized cumulative |dk| evaluated at those indices must reproduce the
    normalized KDIST.  Return the maximum deviation, or None if not evaluable.
    """
    if not labels:
        return None
    n = kpts.shape[0]
    idx = np.array([li.idx_1based - 1 for li in labels], dtype=int)
    if idx.size < 2 or int(idx.min()) < 0 or int(idx.max()) >= n:
        return None

    kds = np.array([li.kdist for li in labels], dtype=float)
    if float(np.nanmax(kds)) == float(np.nanmin(kds)):
        return None

    tgt = normalize_01(kds)
    x = normalize_01(kdist_cumulative_nobreak(kpts, detect_path_breaks(kpts)))
    return float(np.nanmax(np.abs(x[idx] - tgt)))


def choose_kpts_convention(
    kpts: np.ndarray,
    labels: List[LabelInfo],
    B: Optional[np.ndarray],
    tag: str = "EPW",
    debug: bool = False,
) -> Tuple[np.ndarray, str, Optional[float]]:
    """
    Pick between 'k as given' and 'k treated as crystal coordinates and
    converted to cartesian 2 pi/alat', using labelinfo.dat as the referee.
    """
    cands: List[Tuple[str, np.ndarray]] = [("as-given (cartesian 2pi/alat)", kpts)]
    if B is not None:
        cands.append(("crystal -> cartesian 2pi/alat", kpts @ B))

    scored: List[Tuple[str, np.ndarray, Optional[float]]] = []
    for name, kk in cands:
        s = label_index_mismatch(kk, labels)
        scored.append((name, kk, s))
        if debug:
            if s is None:
                print(f"[debug] {tag} k-convention '{name}': not evaluable against labelinfo")
            else:
                print(f"[debug] {tag} k-convention '{name}': max label mismatch = {s:.3e}")

    usable = [t for t in scored if t[2] is not None]
    if not usable:
        return scored[0][1], scored[0][0], None

    name, kk, s = min(usable, key=lambda t: t[2])
    return kk, name, s


def detect_path_breaks(kpts: np.ndarray, factor: float = 3.0) -> List[int]:
    """
    Locate discontinuities in a band k-path, i.e. steps that are far longer
    than the sampling step of the path.  A k-path such as G-X | X-W jumps
    between two symmetry-equivalent images of X; that jump carries no band
    dispersion and must not enter the path length, otherwise everything after
    it is pushed along the normalized axis.

    Returns the 0-based indices of the offending *steps* (between point i and
    point i+1).
    """
    if kpts.shape[0] < 3:
        return []
    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    pos = dk[dk > 0.0]
    if pos.size == 0:
        return []
    med = float(np.median(pos))
    if med <= 0.0:
        return []
    return [int(i) for i in np.where(dk > factor * med)[0]]


def kdist_cumulative_nobreak(kpts: np.ndarray, breaks: Optional[List[int]] = None) -> np.ndarray:
    """Cumulative |dk| with the listed steps counted as zero length."""
    if kpts.shape[0] <= 1:
        return np.zeros(kpts.shape[0], dtype=float)
    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    if breaks:
        dk[np.array(breaks, dtype=int)] = 0.0
    return np.concatenate([[0.0], np.cumsum(dk)])


def labels_to_cart(labels: List[LabelInfo], B: Optional[np.ndarray]) -> Optional[np.ndarray]:
    """High-symmetry points of labelinfo.dat in cartesian 2 pi/alat."""
    if (not labels) or (B is None):
        return None
    if any(li.kfrac is None for li in labels):
        return None
    frac = np.array([li.kfrac for li in labels], dtype=float)
    return frac @ B


def pick_epw_conversion_by_labels(
    kpts: np.ndarray,
    labels: List[LabelInfo],
    lab_cart: np.ndarray,
    B: Optional[np.ndarray],
    tol: float = 1.0e-3,
    debug: bool = False,
) -> Tuple[Optional[np.ndarray], str, Optional[float]]:
    """
    labelinfo.dat gives the index of every high-symmetry point along the band
    k-path, so the k-coordinate convention of band.eig can be read off simply
    by comparing the k-points sitting at those indices with the crystal and
    cartesian coordinates of the labels.
    """
    n = kpts.shape[0]
    idx = np.array([li.idx_1based - 1 for li in labels], dtype=int)
    if idx.size < 2 or int(idx.min()) < 0 or int(idx.max()) >= n:
        if debug:
            print(f"[debug] labelinfo indices out of range for nks={n}")
        return None, "", None

    kk = kpts[idx]
    lab_frac = np.array([li.kfrac for li in labels], dtype=float)
    e_cart = float(np.max(np.abs(kk - lab_cart)))
    e_frac = float(np.max(np.abs(kk - lab_frac)))
    if debug:
        print(f"[debug] EPW k at label indices vs cartesian labels: max dev = {e_cart:.3e}")
        print(f"[debug] EPW k at label indices vs crystal   labels: max dev = {e_frac:.3e}")

    if min(e_cart, e_frac) > tol:
        return None, "", None
    if e_cart <= e_frac:
        return kpts, "as-given (cartesian 2pi/alat)", e_cart
    if B is None:
        return None, "", None
    return kpts @ B, "crystal -> cartesian 2pi/alat", e_frac


def match_labels_in_path(
    kpts: np.ndarray,
    lab_cart: np.ndarray,
    B: Optional[np.ndarray] = None,
    tol: float = 1.0e-4,
    debug: bool = False,
    tag: str = "DFT",
) -> Optional[List[int]]:
    """
    Locate the high-symmetry points inside a k-point list, scanning forward so
    that a repeated point (G at the start and in the middle) is matched in path
    order.  If an exact match fails, points differing by a reciprocal lattice
    vector are accepted as a fallback.
    """
    idxs: List[int] = []
    start = 0
    Binv = np.linalg.inv(B) if B is not None else None

    for i, kl in enumerate(lab_cart):
        sub = kpts[start:]
        if sub.shape[0] == 0:
            return None
        d = np.linalg.norm(sub - kl, axis=1)
        j = int(np.argmin(d))
        if float(d[j]) > tol:
            if Binv is None:
                if debug:
                    print(f"[debug] {tag}: no k-point matches label #{i+1} (min dev {float(d[j]):.3e})")
                return None
            frac = (sub - kl) @ Binv
            r = np.max(np.abs(frac - np.round(frac)), axis=1)
            j = int(np.argmin(r))
            if float(r[j]) > 1.0e-3:
                if debug:
                    print(f"[debug] {tag}: no k-point matches label #{i+1} "
                          f"(min dev {float(d[j]):.3e}, min dev mod G {float(r[j]):.3e})")
                return None
        idxs.append(j + start)
        start = j + start

    return idxs


def report_path_geometry(
    kpts: np.ndarray,
    labels: List[LabelInfo],
    lab_cart: Optional[np.ndarray],
    tag: str = "DFT",
) -> None:
    """
    Debug helper for the case where a dataset cannot be matched to
    labelinfo.dat: print both vertex lists so a genuinely different k-path is
    visible at a glance.
    """
    bounds = detect_path_boundaries(kpts)
    print(f"[debug] {tag} path vertices (cartesian 2pi/alat), {len(bounds)-1} segments:")
    for a, b in zip(bounds, bounds[1:]):
        ka, kb = kpts[a], kpts[b]
        seg = float(np.linalg.norm(kb - ka))
        print(f"[debug]   ({ka[0]:7.4f} {ka[1]:7.4f} {ka[2]:7.4f}) -> "
              f"({kb[0]:7.4f} {kb[1]:7.4f} {kb[2]:7.4f})  |dk| = {seg:.4f}")

    if lab_cart is None or not labels:
        return
    print(f"[debug] labelinfo path vertices (cartesian 2pi/alat), {len(labels)-1} segments:")
    for i in range(len(labels) - 1):
        ka, kb = lab_cart[i], lab_cart[i + 1]
        seg = float(np.linalg.norm(kb - ka))
        dkd = labels[i + 1].kdist - labels[i].kdist
        print(f"[debug]   {labels[i].label:>2s} ({ka[0]:7.4f} {ka[1]:7.4f} {ka[2]:7.4f}) -> "
              f"{labels[i+1].label:>2s} ({kb[0]:7.4f} {kb[1]:7.4f} {kb[2]:7.4f})  "
              f"|dk| = {seg:.4f}, dKDIST = {dkd:.4f}")



def report_corner_positions(
    kpts: np.ndarray,
    labels: List[LabelInfo],
    tag: str = "DFT",
) -> None:
    """
    Debug helper: where the geometric corners of a path sit on the normalized
    x axis, next to the label positions implied by labelinfo KDIST.
    """
    bounds = detect_path_boundaries(kpts)
    x = normalize_01(kdist_cumulative_nobreak(kpts, detect_path_breaks(kpts)))
    xs = [float(x[i]) for i in bounds]
    print(f"[debug] {tag} detected corners (0-based idx -> x): "
          + ", ".join(f"{i}->{v:.6f}" for i, v in zip(bounds, xs)))
    if labels:
        tgt = normalize_01(np.array([li.kdist for li in labels], dtype=float))
        print("[debug] labelinfo positions: "
              + ", ".join(f"{li.label}->{t:.6f}" for li, t in zip(labels, tgt)))
        if len(tgt) == len(xs):
            dev = float(np.nanmax(np.abs(np.array(xs) - tgt)))
            print(f"[debug] {tag} corner-vs-label max deviation = {dev:.3e}")
        else:
            print(f"[debug] {tag} corner count {len(xs)} != label count {len(tgt)} "
                  "(path definitions may differ, or corners are unresolved)")


# ============================================================
# 3) Ef parsing & search
# ============================================================

_FERMI_RE = re.compile(
    r"the\s+Fermi\s+energy\s+is\s+([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*ev",
    re.IGNORECASE,
)

def parse_fermi_from_qe_out(path: str) -> Optional[float]:
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
    for pat in ("*nscf*.out", "*scf*.out"):
        cands = sorted(glob.glob(str(Path(cwd) / pat)))
        if cands:
            return cands[0]
    return None

def resolve_ef(no_fermi_search: bool, fermi_from: Optional[str], set_fermi: Optional[float]) -> Tuple[Optional[float], str]:
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
# 4) PDOS discovery + parsing (QE projwfc outputs)
# ============================================================

_PDOS_ATM_RE = re.compile(
    r"""
    ^(?P<prefix>.+?)\.pdos_atm\#(?P<atm>\d+)\((?P<species>[^)]+)\)_wfc\#(?P<wfc>\d+)\((?P<orb>[^)]+)\)$
    """,
    re.VERBOSE,
)

def _read_numeric_table(path: Path) -> Optional[np.ndarray]:
    """
    Load whitespace-separated numeric table (ignoring comment lines).
    Returns array shape (N, M) or None.
    """
    try:
        rows: List[List[float]] = []
        with open(path, "r") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#"):
                    continue
                parts = s.split()
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    # Skip non-numeric lines silently
                    continue
        if not rows:
            return None
        arr = np.array(rows, dtype=float)
        if arr.ndim != 2 or arr.shape[1] < 2:
            return None
        return arr
    except Exception:
        return None

def find_total_pdos_file(cwd: str = ".") -> Optional[Path]:
    cands = sorted(Path(cwd).glob("*.pdos_tot"))
    return cands[0] if cands else None

@dataclass(frozen=True)
class PdosKey:
    species: str          # e.g. "Ti"
    species_idx: int      # renumbered per species: 1..N
    orbital: str          # e.g. "s", "p", "d", "f"

def discover_projected_pdos_files(prefix: str, cwd: str = ".") -> List[Tuple[Path, int, str, int, str]]:
    """
    Return list of tuples:
      (path, atm_index_original, species, wfc_index, orbital_letter)
    orbital_letter is derived from the '(... )' after wfc#:
      usually 's', 'p', 'd', 'f' (sometimes 'sp' etc). We take the first char.
    """
    out: List[Tuple[Path, int, str, int, str]] = []
    for p in sorted(Path(cwd).glob(f"{prefix}.pdos_atm#*")):
        m = _PDOS_ATM_RE.match(p.name)
        if not m:
            continue
        atm = int(m.group("atm"))
        species = m.group("species").strip()
        wfc = int(m.group("wfc"))
        orb_raw = m.group("orb").strip()
        orb_letter = orb_raw[0].lower() if orb_raw else "?"
        out.append((p, atm, species, wfc, orb_letter))
    return out

def renumber_species_instances(entries: List[Tuple[Path, int, str, int, str]]) -> Dict[Tuple[str, int], int]:
    """
    Build mapping (species, atm_original) -> species_instance_idx (1..)
    For each species, sort by atm_original and assign 1..N.
    """
    by_species: Dict[str, List[int]] = {}
    for _, atm, species, _, _ in entries:
        by_species.setdefault(species, []).append(atm)
    mapping: Dict[Tuple[str, int], int] = {}
    for sp, atms in by_species.items():
        uniq = sorted(set(atms))
        for i, atm in enumerate(uniq, start=1):
            mapping[(sp, atm)] = i
    return mapping

def species_to_cmap_name(species_list: List[str]) -> Dict[str, str]:
    """
    Assign colormaps by species order (stable):
      spring, summer, winter, autumn, spring, ...
    """
    palette = ["spring", "summer", "winter", "autumn"]
    out: Dict[str, str] = {}
    for i, sp in enumerate(species_list):
        out[sp] = palette[i % len(palette)]
    return out

def load_pdos_grouped(
    total_pdos_path: Optional[Path],
    proj_entries: List[Tuple[Path, int, str, int, str]],
    ef: Optional[float],
    align_fermi: bool,
    debug: bool = False,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[PdosKey, np.ndarray]]:
    """
    Returns:
      E_grid (N,), total_dos (N,), proj_dict {PdosKey: pdos(N,)}
    All interpolated onto the total DOS energy grid if possible.
    If total DOS missing, tries to use the first projected file as grid.
    """
    # 1) Decide energy grid
    E_grid: Optional[np.ndarray] = None
    total_dos: Optional[np.ndarray] = None

    if total_pdos_path is not None and total_pdos_path.is_file():
        arr = _read_numeric_table(total_pdos_path)
        if arr is not None:
            E_grid = arr[:, 0].copy()
            # Heuristic: total DOS is usually col 2
            total_dos = arr[:, 1].copy()
            if debug:
                print(f"[debug] total DOS from {total_pdos_path.name}: shape={arr.shape}")

    # Fallback grid from first projected
    if E_grid is None and proj_entries:
        arr0 = _read_numeric_table(proj_entries[0][0])
        if arr0 is not None:
            E_grid = arr0[:, 0].copy()
            if debug:
                print(f"[debug] energy grid from projected file {proj_entries[0][0].name}: shape={arr0.shape}")

    if E_grid is None:
        return None, None, {}

    # Apply fermi shift to energy axis if requested
    if align_fermi and (ef is not None):
        E_grid = E_grid - ef

    # 2) Group projected PDOS by (species_instance, orbital) summing wfcs
    sp_instance_map = renumber_species_instances(proj_entries)
    proj_dict: Dict[PdosKey, np.ndarray] = {}

    for path, atm, species, wfc, orb_letter in proj_entries:
        arr = _read_numeric_table(path)
        if arr is None:
            if debug:
                print(f"[debug] skip unreadable projected PDOS: {path.name}")
            continue

        E = arr[:, 0].copy()
        Y = arr[:, 1:].copy()
        # Combine all PDOS columns for this file (sum m-components etc.)
        y_sum = np.nansum(Y, axis=1)

        # Interpolate onto E_grid
        # First apply Ef shift to E if needed
        if align_fermi and (ef is not None):
            E = E - ef

        # Ensure monotonic for interpolation
        order = np.argsort(E)
        E_sorted = E[order]
        y_sorted = y_sum[order]

        # Remove duplicate E points to avoid np.interp weirdness
        E_unique, idx_unique = np.unique(E_sorted, return_index=True)
        y_unique = y_sorted[idx_unique]

        y_on_grid = np.interp(E_grid, E_unique, y_unique, left=0.0, right=0.0)

        sp_idx = sp_instance_map.get((species, atm), 1)
        key = PdosKey(species=species, species_idx=sp_idx, orbital=orb_letter)

        if key in proj_dict:
            proj_dict[key] = proj_dict[key] + y_on_grid
        else:
            proj_dict[key] = y_on_grid

    # 3) If total_dos exists but on different grid originally, it's already on E_grid
    # If total_dos missing, keep None
    return E_grid, total_dos, proj_dict


# ============================================================
# 5) Plot
# ============================================================

def plot_bands_and_pdos(
    dft_file: str,
    wann_file: str,
    outpng: str,
    ylim: Tuple[float, float],
    ef: Optional[float],
    ef_src: str,
    align_fermi: bool,
    wannier_ef_override: Optional[float],
    labelinfo_file: Optional[str],
    debug: bool = False,
    label_fontsize: int = 16,
    tick_fontsize: int = 14,
    epw_mode: bool = False,
    cell_file: Optional[str] = None,
) -> None:

    # ----- bands -----
    kpts, dft_E = parse_qe_plot_format(dft_file)
    dft_x = normalize_01(kdist_cumulative_raw(kpts))

    epw_tick_xs: Optional[List[float]] = None

    if epw_mode:
        epw_kpts_raw, epw_E = parse_qe_plot_format(wann_file)

        labels_for_x: List[LabelInfo] = []
        if labelinfo_file is not None and Path(labelinfo_file).is_file():
            labels_for_x = read_labelinfo_file(labelinfo_file)

        B = parse_reciprocal_axes(cell_file)
        if debug:
            if B is None:
                print(f"[debug] reciprocal axes b(1..3) not found (cell file: {cell_file}); "
                      "crystal->cartesian conversion unavailable")
            else:
                print(f"[debug] reciprocal axes from {cell_file} (2 pi/alat):")
                for i in range(3):
                    print(f"[debug]   b({i+1}) = "
                          f"({B[i, 0]:.6f} {B[i, 1]:.6f} {B[i, 2]:.6f})")
            print(f"[debug] nks: DFT={kpts.shape[0]}, EPW={epw_kpts_raw.shape[0]}, "
                  f"labels={len(labels_for_x)}")

        lab_cart = labels_to_cart(labels_for_x, B)
        anchored = False

        # Preferred route: labelinfo.dat pins both datasets to the same axis.
        # It carries the index and the physical path length KDIST of every
        # high-symmetry point, so each dataset can be mapped segment by segment
        # onto the normalized KDIST.  This is immune both to a different
        # k-coordinate convention and to path discontinuities (G-X | X-W),
        # whose jump would otherwise be counted as real path length.
        if lab_cart is not None:
            epw_kpts, epw_conv, epw_dev = pick_epw_conversion_by_labels(
                epw_kpts_raw, labels_for_x, lab_cart, B, debug=debug
            )
            epw_idx = [li.idx_1based - 1 for li in labels_for_x]
            dft_idx = match_labels_in_path(kpts, lab_cart, B, debug=debug, tag="DFT")

            ok_epw = (epw_kpts is not None
                      and epw_idx[0] == 0
                      and epw_idx[-1] == epw_kpts.shape[0] - 1
                      and all(b >= a for a, b in zip(epw_idx, epw_idx[1:])))
            ok_dft = (dft_idx is not None
                      and dft_idx[0] == 0
                      and dft_idx[-1] == kpts.shape[0] - 1)

            if debug:
                if epw_kpts is not None:
                    print(f"[debug] EPW k-convention: {epw_conv} (max dev {epw_dev:.3e})")
                print(f"[debug] EPW label indices (0-based): {epw_idx} -> usable={ok_epw}")
                print(f"[debug] DFT label indices (0-based): {dft_idx} -> usable={ok_dft}")

            if ok_epw and ok_dft:
                anchors = normalize_01(np.array([li.kdist for li in labels_for_x], dtype=float))
                dft_x = piecewise_common_x(kpts, list(dft_idx), anchors)
                epw_x = piecewise_common_x(epw_kpts, list(epw_idx), anchors)
                anchored = True
                if debug:
                    print("[debug] x axis anchored on labelinfo KDIST: "
                          + ", ".join(f"{li.label}->{a:.6f}"
                                      for li, a in zip(labels_for_x, anchors)))

        if not anchored:
            # Fallback: no usable labelinfo geometry.  Normalize each path by
            # its own cumulative |dk|, but drop discontinuity jumps first.
            epw_kpts, epw_conv, epw_score = choose_kpts_convention(
                epw_kpts_raw, labels_for_x, B, tag="EPW", debug=debug
            )
            dft_breaks = detect_path_breaks(kpts)
            epw_breaks = detect_path_breaks(epw_kpts)
            if debug:
                if lab_cart is not None:
                    report_path_geometry(kpts, labels_for_x, lab_cart, tag="DFT")
                print(f"[debug] fallback normalization; detected breaks: "
                      f"DFT steps {dft_breaks}, EPW steps {epw_breaks}")
            dft_x = normalize_01(kdist_cumulative_nobreak(kpts, dft_breaks))
            epw_x = normalize_01(kdist_cumulative_nobreak(epw_kpts, epw_breaks))

            if (epw_score is None) or (epw_score > 2.0e-2):
                dft_bounds = detect_path_boundaries(kpts)
                epw_bounds = detect_path_boundaries(epw_kpts)
                if len(dft_bounds) == len(epw_bounds) and len(epw_bounds) >= 2:
                    anchors = normalize_01(kdist_cumulative_nobreak(epw_kpts, epw_breaks))
                    anchors = anchors[np.array(epw_bounds, dtype=int)]
                    dft_x = piecewise_common_x(kpts, dft_bounds, anchors)
                    epw_x = piecewise_common_x(epw_kpts, epw_bounds, anchors)
                    epw_tick_xs = anchors.tolist()
                    if debug:
                        print(f"[debug] DFT path boundaries: {dft_bounds}")
                        print(f"[debug] EPW path boundaries: {epw_bounds}")
                        print(f"[debug] common x anchors: {epw_tick_xs}")
                elif debug:
                    print(
                        "[debug] WARNING: DFT/EPW path segment counts differ: "
                        f"DFT={len(dft_bounds)-1}, EPW={len(epw_bounds)-1}. "
                        "Using independent x normalization; check the k-path input."
                    )
                    report_corner_positions(kpts, labels_for_x, tag="DFT")
                    report_corner_positions(epw_kpts, labels_for_x, tag="EPW")

        wx_blocks = [epw_x.copy() for _ in range(epw_E.shape[0])]
        wy_blocks = [epw_E[ib, :].copy() for ib in range(epw_E.shape[0])]
    else:
        wx_blocks, wy_blocks = parse_wannier_2col_blocks(wann_file)
        wx_blocks = normalize_wannier_x_blocks(wx_blocks)

    if align_fermi and (ef is not None):
        dft_E = dft_E - ef
        ef_w = ef if wannier_ef_override is None else float(wannier_ef_override)
        wy_blocks = [yb - ef_w for yb in wy_blocks]

    # ----- PDOS discovery -----
    total_pdos_path = find_total_pdos_file(".")
    prefix_for_pdos: Optional[str] = None
    if total_pdos_path is not None:
        prefix_for_pdos = total_pdos_path.name.replace(".pdos_tot", "")

    proj_entries: List[Tuple[Path, int, str, int, str]] = []
    if prefix_for_pdos is not None:
        proj_entries = discover_projected_pdos_files(prefix_for_pdos, ".")

    E_dos, total_dos, proj_dict = load_pdos_grouped(
        total_pdos_path=total_pdos_path,
        proj_entries=proj_entries,
        ef=ef,
        align_fermi=align_fermi,
        debug=debug,
    )

    have_dos_panel = (E_dos is not None) and ((total_dos is not None) or (len(proj_dict) > 0))

    # ----- figure layout: 4:1 -----
    fig = plt.figure(figsize=(9, 6))
    if have_dos_panel:
        gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax_band = fig.add_subplot(gs[0, 0])
        ax_dos = fig.add_subplot(gs[0, 1], sharey=ax_band)
    else:
        gs = fig.add_gridspec(1, 1)
        ax_band = fig.add_subplot(gs[0, 0])
        ax_dos = None

    # ----- plot bands (left) -----
    for b in range(dft_E.shape[0]):
        label = "DFT (QE)" if b == 0 else None
        ax_band.plot(dft_x, dft_E[b], color="0.4", linewidth=2.0, alpha=0.6, label=label)

    for i, (xb, yb) in enumerate(zip(wx_blocks, wy_blocks)):
        label = ("EPW interpolation" if epw_mode else "Wannier (MLWF)") if i == 0 else None
        ax_band.plot(xb, yb, "r--", linewidth=1.5, alpha=1.0, label=label)

    # ---- High-symmetry labels ----
    if labelinfo_file is not None:
        p = Path(labelinfo_file)
        if p.is_file():
            labels = read_labelinfo_file(str(p), debug=debug)
            xs, labs = label_x_positions_from_kdist(labels, debug=debug)

            # In EPW mode, if path segmentation was detected successfully,
            # use the common segment anchors for tick positions.  Keep the
            # label strings from labelinfo.dat, but do not use its KDIST to
            # define a second, inconsistent x-axis metric.
            if epw_mode and epw_tick_xs is not None and len(labs) == len(epw_tick_xs):
                xs = epw_tick_xs

            if xs:
                for x in xs:
                    ax_band.axvline(x, color="0.75", linewidth=1.0, linestyle="-", alpha=0.9, zorder=0)
                ax_band.set_xticks(xs)
                ax_band.set_xticklabels(labs, fontsize=tick_fontsize)
        else:
            if debug:
                print(f"[debug] labelinfo file not found: {labelinfo_file}")

    ax_band.set_xlim(0.0, 1.0)
    ax_band.set_ylim(ylim[0], ylim[1])
    ax_band.set_xlabel("Normalized Path", fontsize=label_fontsize)
    ax_band.set_ylabel("Energy (eV)", fontsize=label_fontsize)
    ax_band.grid(True, linestyle=":", alpha=0.6)
    ax_band.tick_params(axis="y", labelsize=tick_fontsize)
    ax_band.legend(loc="upper left")

    # ----- plot DOS (right) -----
    if ax_dos is not None and E_dos is not None:
        # Total DOS in red
        max_x = 0.0
        if total_dos is not None:
            ax_dos.plot(total_dos, E_dos, color="red", linewidth=1.6, alpha=0.95, label="Total DOS")
            max_x = max(max_x, float(np.nanmax(total_dos)))

        # Projected DOS coloring
        # Species order: stable sort by species name as it appears
        species_order = sorted({k.species for k in proj_dict.keys()})
        cmap_by_species = species_to_cmap_name(species_order)

        # For each species-instance, determine orbital list and assign evenly spaced colors
        # Group keys by (species, species_idx)
        by_inst: Dict[Tuple[str, int], List[PdosKey]] = {}
        for key in proj_dict.keys():
            by_inst.setdefault((key.species, key.species_idx), []).append(key)

        # Stable ordering: species then instance idx then orbital
        for (sp, sp_i) in sorted(by_inst.keys(), key=lambda t: (t[0], t[1])):
            keys = by_inst[(sp, sp_i)]
            # Unique orbitals sorted s,p,d,f,... (lexicographic is fine here)
            keys_sorted = sorted(keys, key=lambda k: k.orbital)

            cmap = plt.get_cmap(cmap_by_species.get(sp, "viridis"))
            n = max(1, len(keys_sorted))

            # sample colors away from very ends for readability
            sample_points = np.linspace(0.15, 0.95, n)

            for j, key in enumerate(keys_sorted):
                color = cmap(sample_points[j])
                y = proj_dict[key]
                label = f"{key.species}#{key.species_idx} {key.orbital}"
                ax_dos.plot(y, E_dos, color=color, linewidth=1.2, alpha=0.95, label=label)
                max_x = max(max_x, float(np.nanmax(y)))

        if max_x <= 0.0 or (not np.isfinite(max_x)):
            max_x = 1.0

        ax_dos.set_xlim(0.0, 1.05 * max_x)
        ax_dos.grid(True, linestyle=":", alpha=0.4)
        ax_dos.set_xlabel("DOS / PDOS", fontsize=label_fontsize)

        # Hide duplicated y tick labels on the right
        plt.setp(ax_dos.get_yticklabels(), visible=False)
        ax_dos.tick_params(axis="x", labelsize=tick_fontsize)

        # Legend: keep compact on right
        ax_dos.legend(
            loc="upper right",
            fontsize=max(8, tick_fontsize - 4),
            frameon=True,
            borderpad=0.3,
            handlelength=1.5,
        )

    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    print(f"Saved to {outpng}")
    if debug:
        print(f"[debug] Ef source: {ef_src}, align_fermi={align_fermi}, Ef={ef}")


def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("dft_band_file")
    p.add_argument("wann_band_file")
    p.add_argument("--out", default="band_comparison.png")
    p.add_argument("--ylim", nargs=2, type=float, default=[-3, 3])
    p.add_argument("--no-fermi-search", action="store_true")
    p.add_argument("--fermi-from", default=None)
    p.add_argument("--set-fermi", type=float, default=None)
    p.add_argument("--no-align-fermi", action="store_true")
    p.add_argument("--wannier-fermi", type=float, default=None)
    p.add_argument("--epw", action="store_true",
                   help="Treat the second band file as EPW band.eig (&plot format) instead of Wannier90 2-column data.")
    p.add_argument("--cell-from", default=None,
                   help="QE output used to read the reciprocal axes b(1..3). "
                        "Defaults to --fermi-from, then an auto-detected scf/nscf output.")
    p.add_argument("--labelinfo", default=None)
    p.add_argument("--debug", action="store_true")
    p.add_argument("--label-fontsize", type=int, default=16)
    p.add_argument("--tick-fontsize", type=int, default=14)

    args = p.parse_args(argv)

    # Default labelinfo auto-detect (same behavior as before)
    labelinfo_file = args.labelinfo
    if labelinfo_file is None:
        stem = Path(args.wann_band_file).name
        cand = Path(f"{Path(stem).stem}.labelinfo.dat")
        labelinfo_file = str(cand) if cand.is_file() else None

    ef, ef_src = resolve_ef(args.no_fermi_search, args.fermi_from, args.set_fermi)
    align_fermi = (not args.no_align_fermi)

    cell_file = args.cell_from
    if cell_file is None:
        cell_file = args.fermi_from if args.fermi_from is not None else find_qe_out_for_fermi(".")

    plot_bands_and_pdos(
        dft_file=args.dft_band_file,
        wann_file=args.wann_band_file,
        outpng=args.out,
        ylim=(args.ylim[0], args.ylim[1]),
        ef=ef,
        ef_src=ef_src,
        align_fermi=align_fermi,
        wannier_ef_override=args.wannier_fermi,
        labelinfo_file=labelinfo_file,
        debug=args.debug,
        label_fontsize=args.label_fontsize,
        tick_fontsize=args.tick_fontsize,
        epw_mode=args.epw,
        cell_file=cell_file,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
