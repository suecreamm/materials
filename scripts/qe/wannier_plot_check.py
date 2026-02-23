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

            out.append(LabelInfo(label=lab, idx_1based=idx, kdist=kd))

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
) -> None:

    # ----- bands -----
    kpts, dft_E = parse_qe_plot_format(dft_file)
    dft_x = normalize_01(kdist_cumulative_raw(kpts))

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
        label = "Wannier (MLWF)" if i == 0 else None
        ax_band.plot(xb, yb, "r--", linewidth=1.5, alpha=1.0, label=label)

    # ---- High-symmetry labels ----
    if labelinfo_file is not None:
        p = Path(labelinfo_file)
        if p.is_file():
            labels = read_labelinfo_file(str(p), debug=debug)
            xs, labs = label_x_positions_from_kdist(labels, debug=debug)

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
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
