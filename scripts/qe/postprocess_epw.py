#!/usr/bin/env python3
"""
postprocess_epw.py (v3.4)

EPW postprocessing for QE/EPW outputs:
- Agg backend (headless)
- Save PNG (transparent) + PDF (white)
- omega unit default: meV (optional THz, cm-1)
- jet colormap used consistently
- Meaningful filenames
- decay.* visualizations included (real-space decay diagnostics)
- Larger figures + rcParams + constrained_layout to prevent label overlap
- Eliashberg outputs:
  - PREFIX.imag_iso_*  (w [eV], znorm(w), delta(w) [eV])
  - PREFIX.pade_iso_*  (if numeric table)

Extra (v3.4):
- FS maps enforce equal x/y scaling robustly (aspect='equal', adjustable='box')
- Optional forcing equal x/y limits for unambiguous square scaling
- Optional 1st Brillouin zone (BZ) boundary overlay on FS maps:
  - Reads reciprocal lattice vectors from ./tmp/{prefix}.save/data-file-schema.xml
  - If XML missing or parsing fails: silently skip BZ overlay (plots still generated)
  - k-point convention handling:
      - We attempt a best-effort 'auto' mode:
          if max(|k|) <= 1.2 -> treat as crystal (fractional of b1,b2,b3)
          else treat as cartesian
      - If BZ overlay looks wrong, set KPOINT_MODE below to "crystal" or "cart".

Usage:
  python postprocess_epw.py PREFIX
  python postprocess_epw.py PREFIX --omega-unit THz
  python postprocess_epw.py PREFIX --omega-unit cm-1

Outputs:
  ./plots_open/{PREFIX}__*.png  (transparent)
  ./plots_open/{PREFIX}__*.pdf  (white)
"""

from __future__ import annotations

import argparse
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =========================
# User-tunable knobs
# =========================
# If BZ overlay looks wrong, change "auto" -> "crystal" or "cart".
KPOINT_MODE = "auto"  # {"auto", "crystal", "cart"}

# When aspect_equal is enabled, also force x/y limits to have identical span
# (helps keep the plot visually square even with colorbars/layout).
FORCE_EQUAL_LIMITS = True

# Neighbor range for Wigner–Seitz construction in reciprocal lattice.
# N=2 is usually enough; if BZ looks wrong/degenerate, try N=3.
BZ_NEIGHBOR_RANGE = 2


# =========================
# Global plotting defaults
# =========================
plt.rcParams.update({
    # Canvas / resolution
    "figure.figsize": (8.5, 6.5),
    "figure.dpi": 200,
    "savefig.dpi": 200,

    # Fonts
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,

    # Other
    "axes.grid": False,
})


# =========================
# Units
# =========================

def convert_omega_from_meV(w_mev: np.ndarray, unit: str) -> Tuple[np.ndarray, str]:
    """Convert omega axis from meV to requested unit."""
    if unit == "meV":
        return w_mev, "meV"
    if unit == "THz":
        return w_mev / 4.135667, "THz"
    if unit == "cm-1":
        return w_mev * 8.065544, "cm$^{-1}$"
    raise ValueError(f"Unknown omega unit: {unit}")


# =========================
# Read tables + header sniffing
# =========================

_NUM_RE = re.compile(r"^[\s]*[+-]?(\d+(\.\d*)?|\.\d+)([eEdD][+-]?\d+)?([\s]+|$)")


def _is_numeric_line(s: str) -> bool:
    return bool(_NUM_RE.match(s.strip().replace("D", "E").replace("d", "e")))


def read_table_with_header(path: Path) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Read whitespace numeric table from file.

    Returns
    -------
    data : (N, M) float ndarray
    labels : list[str]
        If a header line provides exactly M tokens, use them; otherwise 'col1..'.
    header_lines : list[str]
        Top comment/header lines (for heuristics).
    """
    lines = path.read_text(errors="replace").splitlines()

    header_lines: List[str] = []
    header_tokens_candidate: Optional[List[str]] = None
    numeric_rows: List[List[float]] = []

    # collect top header/comment lines
    for i in range(min(len(lines), 80)):
        s = lines[i].strip()
        if not s:
            continue
        if _is_numeric_line(s):
            break
        if s.startswith("#") or s.startswith("!") or any(c.isalpha() for c in s):
            header_lines.append(s)
            toks = s.lstrip("#!;@").split()
            if toks and any(any(ch.isalpha() for ch in t) for t in toks):
                header_tokens_candidate = toks

    # parse numeric rows
    for s in lines:
        ss = s.strip()
        if not ss:
            continue
        if ss.startswith(("#", "!", ";", "@")):
            continue
        if not _is_numeric_line(ss):
            continue
        parts = ss.replace("D", "E").replace("d", "e").split()
        try:
            vals = [float(x) for x in parts]
        except ValueError:
            continue
        numeric_rows.append(vals)

    if not numeric_rows:
        raise ValueError(f"No numeric data found in {path}")

    ncol = min(len(r) for r in numeric_rows)
    data = np.array([r[:ncol] for r in numeric_rows], dtype=float)

    if header_tokens_candidate and len(header_tokens_candidate) == ncol:
        labels = header_tokens_candidate
    else:
        labels = [f"col{i}" for i in range(1, ncol + 1)]

    return data, labels, header_lines


# =========================
# Save figures
# =========================

def save_fig(fig: plt.Figure, outbase: Path, dpi: int = 200) -> None:
    """Save PNG (transparent) and PDF (white) with consistent tight bounding."""
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outbase.with_suffix(".png"), dpi=dpi, bbox_inches="tight", transparent=True)
    fig.savefig(outbase.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight",
                transparent=False, facecolor="white")
    plt.close(fig)


def jet_colors(n: int):
    if n <= 1:
        return [plt.cm.jet(0.0)]
    xs = np.linspace(0.0, 1.0, n)
    return [plt.cm.jet(x) for x in xs]


# =========================
# QE XML -> Reciprocal lattice -> 2D 1st BZ polygon (Wigner–Seitz)
# =========================

def _strip_ns(tag: str) -> str:
    """Strip XML namespace from tag."""
    return tag.split("}", 1)[-1] if "}" in tag else tag


def read_reciprocal_lattice_from_qe_xml(xml_path: Path) -> np.ndarray:
    """
    Read reciprocal lattice vectors from QE data-file-schema.xml.

    Returns
    -------
    b : (3, 3) ndarray
        Reciprocal lattice vectors as rows [b1; b2; b3].
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    candidates = {
        "RECIPROCAL_LATTICE_VECTORS",
        "reciprocal_lattice_vectors",
        "RECIPROCAL_LATTICE",
        "reciprocal_lattice",
    }

    node = None
    for elem in root.iter():
        if _strip_ns(elem.tag) in candidates:
            node = elem
            break
    if node is None:
        raise ValueError("Reciprocal lattice tag not found in XML")

    vecs = []
    for child in list(node):
        txt = (child.text or "").strip()
        if not txt:
            continue
        parts = txt.replace("D", "E").replace("d", "e").split()
        if len(parts) >= 3:
            vecs.append([float(parts[0]), float(parts[1]), float(parts[2])])

    # Fallback: try to pull 9 floats from within node
    if len(vecs) != 3:
        all_txt = " ".join([(child.text or "") for child in list(node)] + [(node.text or "")])
        toks = all_txt.replace("D", "E").replace("d", "e").split()
        vals = []
        for t in toks:
            try:
                vals.append(float(t))
            except ValueError:
                pass
        if len(vals) >= 9:
            vecs = [vals[0:3], vals[3:6], vals[6:9]]

    if len(vecs) != 3:
        raise ValueError(f"Could not parse 3 reciprocal vectors (parsed={len(vecs)})")

    return np.array(vecs, dtype=float)


def _clip_polygon_halfplane(poly: np.ndarray, n: np.ndarray, c: float, eps: float = 1e-12) -> np.ndarray:
    """
    Clip a convex polygon by half-plane: n·x <= c

    Parameters
    ----------
    poly : (M,2) array
        Polygon vertices in CCW order.
    n : (2,) array
        Half-plane normal.
    c : float
        Half-plane offset.

    Returns
    -------
    new_poly : (K,2) array
        Clipped polygon vertices. May be empty.
    """
    if poly.size == 0:
        return poly

    def inside(p):
        return (n @ p) <= (c + eps)

    def intersect(p1, p2):
        d = p2 - p1
        denom = n @ d
        if abs(denom) < eps:
            return p1
        t = (c - (n @ p1)) / denom
        t = max(0.0, min(1.0, float(t)))
        return p1 + t * d

    out = []
    prev = poly[-1]
    prev_in = inside(prev)

    for cur in poly:
        cur_in = inside(cur)
        if cur_in:
            if not prev_in:
                out.append(intersect(prev, cur))
            out.append(cur)
        else:
            if prev_in:
                out.append(intersect(prev, cur))
        prev, prev_in = cur, cur_in

    return np.array(out, dtype=float)


def first_bz_polygon_2d_from_b1b2(b1: np.ndarray, b2: np.ndarray, N: int = 2) -> np.ndarray:
    """
    Compute 2D 1st Brillouin zone as Wigner–Seitz cell of reciprocal lattice generated by b1,b2.

    Method:
    - Start with a large bounding square.
    - For each nearby reciprocal lattice point G != 0, clip by the half-plane:
        x · G <= |G|^2 / 2
      which keeps points closer to the origin than to G.

    Returns
    -------
    poly : (M,2) array
        Vertices of the 1st BZ polygon in CCW order.
    """
    Gs = []
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if m == 0 and n == 0:
                continue
            G = m * b1 + n * b2
            if np.linalg.norm(G) > 0:
                Gs.append(G)

    if not Gs:
        raise ValueError("No neighbor G vectors generated")

    R = 2.5 * max(np.linalg.norm(g) for g in Gs)
    poly = np.array([[-R, -R], [R, -R], [R, R], [-R, R]], dtype=float)

    for G in Gs:
        nvec = G
        c = 0.5 * float(G @ G)
        poly = _clip_polygon_halfplane(poly, nvec, c)
        if poly.size == 0:
            break

    if poly.size == 0 or poly.shape[0] < 3:
        raise ValueError("BZ polygon clipping produced empty/degenerate polygon")

    return poly


def kpoints_to_cart_2d(kxyz: np.ndarray, b: np.ndarray, mode: str = "auto") -> np.ndarray:
    """
    Convert k-points to 2D Cartesian for plotting/overlay consistency.

    Parameters
    ----------
    kxyz : (N,3) array
        k-point columns as read from lambda_FS.
    b : (3,3) array
        Reciprocal lattice vectors as rows [b1;b2;b3].
    mode : {"auto","crystal","cart"}
        - "crystal": interpret kxyz as fractional coefficients along (b1,b2,b3),
          i.e., k_cart = k1*b1 + k2*b2 + k3*b3.
        - "cart": interpret kxyz already as Cartesian in the same basis as b (no transform).
        - "auto": heuristic:
            if max(|k|) <= 1.2 -> treat as crystal; else treat as cart.

    Returns
    -------
    kxy : (N,2) array
        2D Cartesian k-points.
    """
    if kxyz.ndim != 2 or kxyz.shape[1] < 3:
        raise ValueError("kxyz must be (N,3)")

    mode_use = mode
    if mode == "auto":
        m = float(np.nanmax(np.abs(kxyz[:, :3])))
        mode_use = "crystal" if (m <= 1.2) else "cart"

    if mode_use == "cart":
        return kxyz[:, :2].copy()

    if mode_use == "crystal":
        k_cart = kxyz[:, :3] @ b  # (N,3)
        return k_cart[:, :2]

    raise ValueError(f"Unknown k-point mode: {mode}")


# =========================
# Plot primitives
# =========================

def plot_multi_line(x: np.ndarray, ys: List[np.ndarray], curve_labels: Optional[List[str]],
                    title: str, xlabel: str, ylabel: str,
                    outbase: Path,
                    figsize: Optional[Tuple[float, float]] = None) -> None:
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    cols = jet_colors(len(ys))
    for i, y in enumerate(ys):
        lab = curve_labels[i] if curve_labels and i < len(curve_labels) else None
        ax.plot(x, y, color=cols[i], label=lab)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if curve_labels and len(ys) > 1:
        ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    save_fig(fig, outbase)


def plot_scatter_2d(x: np.ndarray, y: np.ndarray, v: np.ndarray,
                    title: str, xlabel: str, ylabel: str, cbar_label: str,
                    outbase: Path,
                    figsize: Tuple[float, float] = (8.2, 6.8),
                    aspect_equal: bool = True,
                    force_equal_limits: bool = True,
                    bz_poly: Optional[np.ndarray] = None) -> None:
    """
    Scatter plot with colormap and optional BZ boundary overlay.

    Parameters
    ----------
    bz_poly:
        (M,2) vertices of BZ polygon in the SAME coordinate system as x,y.
        If provided, it will be drawn as a closed line.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    sc = ax.scatter(x, y, c=v, s=10, cmap="jet")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if aspect_equal:
        ax.set_aspect("equal", adjustable="box")

        if force_equal_limits and (x.size > 0) and (y.size > 0):
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            ymin, ymax = np.nanmin(y), np.nanmax(y)
            xmid, ymid = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
            span = max(xmax - xmin, ymax - ymin)
            pad = 0.02 * span if span > 0 else 0.0
            ax.set_xlim(xmid - 0.5 * span - pad, xmid + 0.5 * span + pad)
            ax.set_ylim(ymid - 0.5 * span - pad, ymid + 0.5 * span + pad)

    # ---- BZ overlay ----
    if bz_poly is not None and isinstance(bz_poly, np.ndarray) and bz_poly.ndim == 2 and bz_poly.shape[1] == 2:
        px = np.r_[bz_poly[:, 0], bz_poly[0, 0]]
        py = np.r_[bz_poly[:, 1], bz_poly[0, 1]]
        ax.plot(px, py, lw=1.2)

    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(cbar_label)
    save_fig(fig, outbase)


def plot_xy_line_plus_points(x: np.ndarray, y: np.ndarray,
                             title: str, xlabel: str, ylabel: str,
                             outbase: Path,
                             figsize: Optional[Tuple[float, float]] = None) -> None:
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.plot(x, y, color=plt.cm.jet(0.2), lw=1.6)
    ax.scatter(x, y, s=14, color=plt.cm.jet(0.7))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    save_fig(fig, outbase)


# =========================
# File-specific plotters
# =========================

def plot_a2f(prefix: str, path: Path, outdir: Path, omega_unit: str) -> None:
    data, labels, _ = read_table_with_header(path)
    if data.shape[1] < 2:
        raise ValueError("a2f needs >=2 columns")
    w_mev = data[:, 0]
    w, wlab = convert_omega_from_meV(w_mev, omega_unit)

    plot_multi_line(
        w, [data[:, 1]],
        curve_labels=[labels[1]] if labels else None,
        title=f"{prefix}: α²F(ω)",
        xlabel=f"ω ({wlab})",
        ylabel="α²F(ω)",
        outbase=outdir / f"{prefix}__a2f__alpha2F_vs_omega",
        figsize=(9.0, 6.5),
    )

    if data.shape[1] >= 3:
        plot_multi_line(
            w, [data[:, 2]],
            curve_labels=[labels[2]] if labels else None,
            title=f"{prefix}: λ(ω) (from a2f col3)",
            xlabel=f"ω ({wlab})",
            ylabel="λ(ω)",
            outbase=outdir / f"{prefix}__a2f__lambda_vs_omega",
            figsize=(9.0, 6.5),
        )


def plot_spectrum(prefix: str, path: Path, outdir: Path, omega_unit: str,
                  tag: str, ylabel: str, title: str, outname: str) -> None:
    data, labels, _ = read_table_with_header(path)
    if data.shape[1] < 2:
        raise ValueError(f"{tag} needs >=2 columns")

    w_mev = data[:, 0]
    w, wlab = convert_omega_from_meV(w_mev, omega_unit)

    ys = [data[:, i] for i in range(1, data.shape[1])]
    curve_labels = [labels[i] for i in range(1, data.shape[1])] if labels else None

    plot_multi_line(
        w, ys,
        curve_labels=curve_labels if (curve_labels and len(ys) > 1) else None,
        title=title,
        xlabel=f"ω ({wlab})",
        ylabel=ylabel,
        outbase=outdir / f"{prefix}__{outname}",
        figsize=(9.0, 6.5),
    )


def plot_lambda_fs(prefix: str, path: Path, outdir: Path, tag: str) -> None:
    """
    Typical formats include:
      kx ky kz band (Enk-Ef) lambda

    We plot:
      - FS map of lambda: color=λ at (kx,ky)
      - optionally Enk-Ef map if present
      - optionally overlay 1st BZ boundary if QE XML is available
    """
    data, _, _ = read_table_with_header(path)
    if data.shape[1] < 3:
        raise ValueError(f"{tag} needs >=3 columns")

    # k-point in file is typically 3 columns
    kxyz = data[:, 0:3]
    lam = data[:, -1]

    # ---- BZ overlay attempt (XML) ----
    bz_poly = None
    kxy = None
    xml_path = Path("./tmp") / f"{prefix}.save" / "data-file-schema.xml"

    if xml_path.exists():
        try:
            b = read_reciprocal_lattice_from_qe_xml(xml_path)
            b1_2d = b[0, :2]
            b2_2d = b[1, :2]
            bz_poly = first_bz_polygon_2d_from_b1b2(b1_2d, b2_2d, N=BZ_NEIGHBOR_RANGE)
            kxy = kpoints_to_cart_2d(kxyz, b, mode=KPOINT_MODE)
            print(f"[info] BZ overlay enabled ({xml_path}) | kpoint_mode={KPOINT_MODE}")
        except Exception as e:
            print(f"[skip] BZ overlay: failed from {xml_path}: {e}")
            bz_poly = None
            kxy = None
    else:
        print(f"[skip] BZ overlay: XML not found at {xml_path}")

    # If we have kxy, use it; otherwise raw columns 0/1
    if kxy is None:
        kx = data[:, 0]
        ky = data[:, 1]
    else:
        kx = kxy[:, 0]
        ky = kxy[:, 1]

    # ---- lambda map ----
    mL = np.isfinite(kx) & np.isfinite(ky) & np.isfinite(lam)
    kxL, kyL, lamL = kx[mL], ky[mL], lam[mL]

    plot_scatter_2d(
        kxL, kyL, lamL,
        title=f"{prefix}: {tag} FS map (λ)",
        xlabel="k_x",
        ylabel="k_y",
        cbar_label="λ",
        outbase=outdir / f"{prefix}__{tag}__FSmap_lambda",
        figsize=(8.6, 7.0),
        aspect_equal=True,
        force_equal_limits=FORCE_EQUAL_LIMITS,
        bz_poly=bz_poly,
    )

    # ---- Enk-Ef map ----
    if data.shape[1] >= 5:
        enk_ef = data[:, -2]
        mE = np.isfinite(kx) & np.isfinite(ky) & np.isfinite(enk_ef)
        plot_scatter_2d(
            kx[mE], ky[mE], enk_ef[mE],
            title=f"{prefix}: {tag} FS map (Enk-Ef)",
            xlabel="k_x",
            ylabel="k_y",
            cbar_label="Enk-Ef (eV)",
            outbase=outdir / f"{prefix}__{tag}__FSmap_Enk_minus_Ef",
            figsize=(8.6, 7.0),
            aspect_equal=True,
            force_equal_limits=FORCE_EQUAL_LIMITS,
            bz_poly=bz_poly,
        )


def plot_lambda_distribution(prefix: str, path: Path, outdir: Path) -> None:
    """
    For PREFIX.lambda_k_pairs:
      col1: lambda_nk
      col2: dist_scaled
      col3: dist_unscaled (optional)
    """
    data, _, _ = read_table_with_header(path)
    if data.shape[1] < 2:
        raise ValueError("lambda_k_pairs distribution needs >=2 columns")

    lam = data[:, 0]
    dist_scaled = data[:, 1]

    ys = [dist_scaled]
    labs = ["dist_scaled"]

    if data.shape[1] >= 3:
        ys.append(data[:, 2])
        labs.append("dist_unscaled")

    plot_multi_line(
        lam, ys, curve_labels=labs,
        title=f"{prefix}: ρ(λ_nk) distribution",
        xlabel="λ_nk",
        ylabel="distribution",
        outbase=outdir / f"{prefix}__lambda_k_pairs__rho_lambda",
        figsize=(9.0, 6.5),
    )


def plot_decay_file(prefix: str, path: Path, outdir: Path) -> None:
    """
    decay.* are real-space decay diagnostics.
    Typical format:
      R [Ang]   value
    """
    data, labels, header = read_table_with_header(path)
    if data.shape[1] < 2:
        raise ValueError("decay file needs >=2 columns")

    x = data[:, 0]
    y = data[:, 1]

    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]

    htxt = " ".join(header)
    xlabel = "R (Angstrom)" if ("Ang" in htxt or "Ang]" in htxt or "Angstrom" in htxt) else (labels[0] if labels else "x")

    if "|g" in htxt or ("g(" in htxt and "Ry" in htxt):
        ylabel = "max |g(R)| (Ry)"
        outmetric = "gmax_vs_R"
    elif "dyn" in path.name.lower():
        ylabel = labels[1] if labels else "dynmat metric"
        outmetric = "dynmat_metric_vs_R"
    elif "H" in path.name:
        ylabel = labels[1] if labels else "H metric"
        outmetric = "H_metric_vs_R"
    else:
        ylabel = labels[1] if labels else "value"
        outmetric = "value_vs_R"

    stem = path.name.replace(".", "_")
    outbase = outdir / f"{prefix}__{stem}__{outmetric}"

    plot_xy_line_plus_points(
        x, y,
        title=f"{prefix}: {path.name}",
        xlabel=xlabel,
        ylabel=ylabel,
        outbase=outbase,
        figsize=(8.8, 6.4),
    )


# =========================
# Eliashberg: imag_iso_* and pade_iso_*
# =========================

def _parse_T_from_suffix(name: str, prefix: str, tag: str) -> str:
    """
    Extract temperature-ish suffix from filenames like:
      PREFIX.imag_iso_005.00
      PREFIX.pade_iso_032.86
    Returns '005.00' etc (string). Best-effort.
    """
    base = name
    key = f"{prefix}.{tag}_"
    if key in base:
        return base.split(key, 1)[1]
    if "_" in base:
        return base.split("_")[-1]
    return "unknownT"


def plot_eliashberg_iso(prefix: str, path: Path, outdir: Path) -> None:
    """
    Plot isotropic Eliashberg results from PREFIX.imag_iso_*.

    Expected columns:
      w [eV]   znorm(w)   delta(w) [eV]
    """
    data, _, _ = read_table_with_header(path)
    if data.shape[1] < 3:
        raise ValueError("imag_iso needs >=3 columns: w, znorm, delta")

    w = data[:, 0]  # eV
    z = data[:, 1]
    d = data[:, 2]  # eV

    mw = np.isfinite(w)
    w = w[mw]; z = z[mw]; d = d[mw]

    Ttag = _parse_T_from_suffix(path.name, prefix, "imag_iso")

    mz = np.isfinite(z)
    if np.any(mz):
        plot_multi_line(
            w[mz], [z[mz]], curve_labels=None,
            title=f"{prefix}: Eliashberg isotropic znorm(ω)  (T={Ttag})",
            xlabel="ω (eV)",
            ylabel="znorm(ω)",
            outbase=outdir / f"{prefix}__eliashberg__imag_iso_{Ttag}__znorm_vs_w_eV",
            figsize=(9.2, 6.6),
        )
    else:
        raise ValueError(f"All znorm are NaN/inf in {path.name}")

    md = np.isfinite(d)
    if np.any(md):
        d_meV = d[md] * 1000.0
        plot_multi_line(
            w[md], [d_meV], curve_labels=None,
            title=f"{prefix}: Eliashberg isotropic Δ(ω)  (T={Ttag})",
            xlabel="ω (eV)",
            ylabel="Δ(ω) (meV)",
            outbase=outdir / f"{prefix}__eliashberg__imag_iso_{Ttag}__delta_vs_w_meV",
            figsize=(9.2, 6.6),
        )
    else:
        raise ValueError(f"All delta are NaN/inf in {path.name}")


def plot_pade_iso(prefix: str, path: Path, outdir: Path) -> None:
    """
    Plot PREFIX.pade_iso_* if it is a numeric table.

    Best-effort:
      - first col assumed ω (eV)
      - remaining columns plotted as separate curves
    """
    data, labels, _ = read_table_with_header(path)
    if data.shape[1] < 2:
        raise ValueError("pade_iso needs >=2 columns")

    w = data[:, 0]  # likely eV
    Ttag = _parse_T_from_suffix(path.name, prefix, "pade_iso")

    ys = [data[:, i] for i in range(1, data.shape[1])]
    mw = np.isfinite(w)
    w = w[mw]
    ys = [y[mw] for y in ys]

    ys_clean = []
    for y in ys:
        yy = y.copy()
        yy[~np.isfinite(yy)] = np.nan
        ys_clean.append(yy)

    curve_labels = None
    if labels and len(labels) == data.shape[1]:
        curve_labels = labels[1:]

    plot_multi_line(
        w, ys_clean,
        curve_labels=curve_labels if (curve_labels and len(ys_clean) > 1) else None,
        title=f"{prefix}: Eliashberg pade_iso  (T={Ttag})",
        xlabel="ω (eV)",
        ylabel="value (arb.)",
        outbase=outdir / f"{prefix}__eliashberg__pade_iso_{Ttag}__cols_vs_w",
        figsize=(9.2, 6.6),
    )


# =========================
# Main
# =========================

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--omega-unit", default="meV", choices=["meV", "THz", "cm-1"])
    args = ap.parse_args()

    prefix = args.prefix
    omega_unit = args.omega_unit
    cwd = Path.cwd()

    outdir = cwd / "plots_open"
    outdir.mkdir(exist_ok=True)

    print(f"[info] prefix={prefix}")
    print(f"[info] omega-unit={omega_unit} (default meV if not specified)")
    print(f"[info] output dir={outdir}")
    print(f"[info] FS aspect equal: True | FORCE_EQUAL_LIMITS={FORCE_EQUAL_LIMITS}")
    print(f"[info] BZ overlay: xml=./tmp/{prefix}.save/data-file-schema.xml | KPOINT_MODE={KPOINT_MODE} | BZ_NEIGHBOR_RANGE={BZ_NEIGHBOR_RANGE}")

    # ---- a2f ----
    p = cwd / f"{prefix}.a2f"
    if p.exists():
        try:
            plot_a2f(prefix, p, outdir, omega_unit)
            print("[ok] a2f (+ lambda(omega) if present)")
        except Exception as e:
            print(f"[skip] {p.name}: {e}")

    # ---- spectra ----
    specs = [
        ("a2f_proj", "α²F(ω) projected (arb.)", f"{prefix}: projected α²F(ω)", "a2f_proj__alpha2Fproj_vs_omega"),
        ("phdos", "Phonon DOS (arb.)", f"{prefix}: phonon DOS", "phdos__dos_vs_omega"),
        ("phdos_proj", "Phonon DOS (arb.)", f"{prefix}: phonon DOS (proj)", "phdos_proj__dosproj_vs_omega"),
    ]
    for tag, ylabel, title, outname in specs:
        p = cwd / f"{prefix}.{tag}"
        if p.exists():
            try:
                plot_spectrum(prefix, p, outdir, omega_unit, tag, ylabel, title, outname)
                print(f"[ok] {tag}")
            except Exception as e:
                print(f"[skip] {p.name}: {e}")

    # ---- lambda FS-like ----
    for tag in ["lambda_FS", "lambda_aniso", "lambda_pairs", "lambda.frmsf"]:
        p = cwd / f"{prefix}.{tag}"
        if p.exists():
            try:
                plot_lambda_fs(prefix, p, outdir, tag)
                print(f"[ok] {tag} -> FS maps")
            except Exception as e:
                print(f"[skip] {p.name}: {e}")

    # ---- lambda_k_pairs distribution ----
    p = cwd / f"{prefix}.lambda_k_pairs"
    if p.exists():
        try:
            plot_lambda_distribution(prefix, p, outdir)
            print("[ok] lambda_k_pairs -> rho(lambda) distribution")
        except Exception as e:
            print(f"[skip] {p.name}: {e}")

    # ---- decay.* ----
    for df in sorted(cwd.glob("decay.*")):
        try:
            plot_decay_file(prefix, df, outdir)
            print(f"[ok] {df.name} -> decay plot")
        except Exception as e:
            print(f"[skip] {df.name}: {e}")

    # ---- Eliashberg: imag_iso_* ----
    for f in sorted(cwd.glob(f"{prefix}.imag_iso_*")):
        try:
            plot_eliashberg_iso(prefix, f, outdir)
            print(f"[ok] {f.name} -> Eliashberg imag_iso plots")
        except Exception as e:
            print(f"[skip] {f.name}: {e}")

    # ---- Eliashberg: pade_iso_* ----
    for f in sorted(cwd.glob(f"{prefix}.pade_iso_*")):
        try:
            plot_pade_iso(prefix, f, outdir)
            print(f"[ok] {f.name} -> Eliashberg pade_iso plots")
        except Exception as e:
            print(f"[skip] {f.name}: {e}")

    print("[done] postprocess complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

