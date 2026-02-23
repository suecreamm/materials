#!/usr/bin/env python3
"""
pdos_components_plot.py

Overview plot:
- Total DOS: light red fill, no line
- Grouped PDOS: species-colormap rule (spring/summer/winter/autumn cycling)

Component plots (per PDOS file):
- pdos_sum: light gray fill, no line
- pdos components: SOLID lines only (no markers), distinct colors (tab10)
- Component labels:
    * If component count is 1/3/5/7, assign conventional real-orbital labels
      via a dictionary (editable in one place).
    * Otherwise fallback to pdos1, pdos2, ...

NOTE (important):
- QE projwfc's internal ordering/convention for p/d/f components may differ from the
  conventional cubic-harmonics ordering used in the dictionary below.
  This script labels components by a chosen convention; if you find a mismatch,
  adjust the label order in ORBITAL_LABELS below (single edit point).
"""

from __future__ import annotations

import argparse
import glob
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# -----------------------------
# Fermi parsing (QE outputs)
# -----------------------------
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


# -----------------------------
# PDOS discovery + parsing
# -----------------------------
_PDOS_ATM_RE = re.compile(
    r"^(?P<prefix>.+?)\.pdos_atm\#(?P<atm>\d+)\((?P<species>[^)]+)\)_wfc\#(?P<wfc>\d+)\((?P<orb>[^)]+)\)$"
)

@dataclass(frozen=True)
class PdosFile:
    path: Path
    atm: int
    species: str
    wfc: int
    orb_raw: str
    orb_letter: str

@dataclass(frozen=True)
class PdosKey:
    species: str
    inst: int
    orbital: str


def _read_numeric_table(path: Path) -> Optional[np.ndarray]:
    """Load whitespace-separated numeric table (ignoring # comment lines)."""
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
                    continue
        if not rows:
            return None
        arr = np.array(rows, dtype=float)
        if arr.ndim != 2 or arr.shape[1] < 2:
            return None
        return arr
    except Exception:
        return None


def find_prefix_from_pdos_tot(cwd: str = ".") -> Optional[str]:
    cands = sorted(Path(cwd).glob("*.pdos_tot"))
    if not cands:
        return None
    return cands[0].name.replace(".pdos_tot", "")

def find_total_dos_file(cwd: str = ".") -> Optional[Path]:
    cands = sorted(Path(cwd).glob("*.pdos_tot"))
    return cands[0] if cands else None

def discover_all_projected_pdos(prefix: str, cwd: str = ".") -> List[PdosFile]:
    out: List[PdosFile] = []
    for p in sorted(Path(cwd).glob(f"{prefix}.pdos_atm#*")):
        m = _PDOS_ATM_RE.match(p.name)
        if not m:
            continue
        atm = int(m.group("atm"))
        species = m.group("species").strip()
        wfc = int(m.group("wfc"))
        orb_raw = m.group("orb").strip()
        orb_letter = orb_raw[0].lower() if orb_raw else "?"
        out.append(PdosFile(path=p, atm=atm, species=species, wfc=wfc, orb_raw=orb_raw, orb_letter=orb_letter))
    return out

def renumber_species_instances(pdos_files: List[PdosFile]) -> Dict[Tuple[str, int], int]:
    by_species: Dict[str, List[int]] = {}
    for pf in pdos_files:
        by_species.setdefault(pf.species, []).append(pf.atm)

    mapping: Dict[Tuple[str, int], int] = {}
    for sp, atms in by_species.items():
        uniq = sorted(set(atms))
        for i, atm in enumerate(uniq, start=1):
            mapping[(sp, atm)] = i
    return mapping

def species_instance_counts(pdos_files: List[PdosFile]) -> Dict[str, int]:
    by_species: Dict[str, set] = {}
    for pf in pdos_files:
        by_species.setdefault(pf.species, set()).add(pf.atm)
    return {sp: len(s) for sp, s in by_species.items()}

def species_to_cmap_name(species_list: List[str]) -> Dict[str, str]:
    palette = ["spring", "summer", "winter", "autumn"]
    return {sp: palette[i % len(palette)] for i, sp in enumerate(species_list)}

def format_species_label(species: str, inst: int, counts: Dict[str, int]) -> str:
    n = counts.get(species, 1)
    return species if n <= 1 else f"{species}#{inst}"

def compute_ewin(emin_arg: Optional[float], emax_arg: Optional[float], default_emin: float, default_emax: float) -> Tuple[float, float]:
    emin = default_emin if emin_arg is None else float(emin_arg)
    emax = default_emax if emax_arg is None else float(emax_arg)
    if emax < emin:
        emin, emax = emax, emin
    return emin, emax


def load_grouped_pdos_on_grid(
    total_dos_path: Optional[Path],
    pdos_files: List[PdosFile],
    ef: Optional[float],
    align_fermi: bool,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[PdosKey, np.ndarray], Dict[str, int]]:
    E_grid: Optional[np.ndarray] = None
    total_dos: Optional[np.ndarray] = None

    if total_dos_path is not None and total_dos_path.is_file():
        arr = _read_numeric_table(total_dos_path)
        if arr is not None:
            E_grid = arr[:, 0].copy()
            total_dos = arr[:, 1].copy()

    if E_grid is None and pdos_files:
        arr0 = _read_numeric_table(pdos_files[0].path)
        if arr0 is not None:
            E_grid = arr0[:, 0].copy()

    if E_grid is None:
        return None, None, {}, {}

    if align_fermi and (ef is not None):
        E_grid = E_grid - ef

    inst_map = renumber_species_instances(pdos_files)
    sp_counts = species_instance_counts(pdos_files)

    grouped: Dict[PdosKey, np.ndarray] = {}
    for pf in pdos_files:
        arr = _read_numeric_table(pf.path)
        if arr is None:
            continue

        E = arr[:, 0].copy()
        y_sum = np.nansum(arr[:, 1:], axis=1)  # ldos + pdos comps

        if align_fermi and (ef is not None):
            E = E - ef

        order = np.argsort(E)
        E_sorted = E[order]
        y_sorted = y_sum[order]
        E_unique, idx_unique = np.unique(E_sorted, return_index=True)
        y_unique = y_sorted[idx_unique]

        y_on = np.interp(E_grid, E_unique, y_unique, left=0.0, right=0.0)

        inst = inst_map.get((pf.species, pf.atm), 1)
        key = PdosKey(species=pf.species, inst=inst, orbital=pf.orb_letter)
        grouped[key] = grouped.get(key, 0.0) + y_on

    return E_grid, total_dos, grouped, sp_counts


def build_color_map_for_keys(
    grouped_keys: List[PdosKey],
    sp_counts: Dict[str, int],
) -> Dict[PdosKey, Tuple[float, float, float, float]]:
    species_order = sorted({k.species for k in grouped_keys})
    cmap_by_species = species_to_cmap_name(species_order)

    by_inst: Dict[Tuple[str, int], List[PdosKey]] = {}
    for k in grouped_keys:
        by_inst.setdefault((k.species, k.inst), []).append(k)

    out: Dict[PdosKey, Tuple[float, float, float, float]] = {}
    for (sp, inst), keys in sorted(by_inst.items(), key=lambda t: (t[0][0], t[0][1])):
        keys_sorted = sorted(keys, key=lambda kk: kk.orbital)
        cmap = plt.get_cmap(cmap_by_species.get(sp, "viridis"))
        n = max(1, len(keys_sorted))
        sample_points = np.linspace(0.15, 0.95, n)
        for j, kk in enumerate(keys_sorted):
            out[kk] = cmap(sample_points[j])
    return out


# -----------------------------
# Conventional real-orbital labels (EDIT HERE if order mismatch)
# -----------------------------

ORBITAL_LABELS = {
    # l=0
    1: ["s"],

    # l=1 (p): m=0, 1cos, 1sin  ->  z, x, y
    3: ["p_z", "p_x", "p_y"],

    # l=2 (d): m=0, 1cos, 1sin, 2cos, 2sin
    #  ->  (3z^2-r^2), xz, yz, (x^2-y^2), xy
    5: ["d_{z^2}", "d_{xz}", "d_{yz}", "d_{x^2-y^2}", "d_{xy}"],

    # l=3 (f): m=0, 1cos, 1sin, 2cos, 2sin, 3cos, 3sin
    # Typical polynomial forms:
    #   m=0:            z(5z^2-3r^2)
    #   m=1 cos/sin:    x(5z^2-r^2), y(5z^2-r^2)
    #   m=2 cos/sin:    z(x^2-y^2),  xyz
    #   m=3 cos/sin:    x(x^2-3y^2), y(3x^2-y^2)
    7: [
        "f_{z(5z^2-3r^2)}",
        "f_{x(5z^2-r^2)}",
        "f_{y(5z^2-r^2)}",
        "f_{z(x^2-y^2)}",
        "f_{xyz}",
        "f_{x(x^2-3y^2)}",
        "f_{y(3x^2-y^2)}",
    ],

    # l=4 (g): m=0, 1cos, 1sin, 2cos, 2sin, 3cos, 3sin, 4cos, 4sin
    # Typical polynomial forms:
    #   m=0:            35z^4 - 30z^2 r^2 + 3 r^4
    #   m=1 cos/sin:    x z(7z^2-3r^2), y z(7z^2-3r^2)
    #   m=2 cos/sin:    (x^2-y^2)(7z^2-r^2),  xy(7z^2-r^2)
    #   m=3 cos/sin:    xz(x^2-3y^2), yz(3x^2-y^2)
    #   m=4 cos/sin:    (x^4-6x^2y^2+y^4),  xy(x^2-y^2)
    9: [
        "g_{35z^4-30z^2r^2+3r^4}",
        "g_{xz(7z^2-3r^2)}",
        "g_{yz(7z^2-3r^2)}",
        "g_{(x^2-y^2)(7z^2-r^2)}",
        "g_{xy(7z^2-r^2)}",
        "g_{xz(x^2-3y^2)}",
        "g_{yz(3x^2-y^2)}",
        "g_{x^4-6x^2y^2+y^4}",
        "g_{xy(x^2-y^2)}",
    ],
}


def component_labels(ncomp: int, orb_letter: str) -> list[str]:
    """
    Return component labels for PDOS columns based on QE ylmr2 ordering.

    Notes
    -----
    - Ordering is fixed by QE `upflib/ylmr2.f90`:
      m=0, (m=1 cos/sin), (m=2 cos/sin), ...
    - Labels are conventional polynomial-shaped names for real (tesseral)
      spherical harmonics. Overall normalization/sign does not matter for
      PDOS component identification.
    """
    base = ORBITAL_LABELS.get(ncomp)
    if base is None:
        return [f"pdos{i}" for i in range(1, ncomp + 1)]

    ol = orb_letter.lower()
    if ol in ("s", "p", "d", "f", "g"):
        # optional prefix to remind which l-shell file it came from
        return [f"{ol}:{lab}" for lab in base]
    return base

# -----------------------------
# Plotting
# -----------------------------
def plot_overview(
    outpath: Path,
    E: np.ndarray,
    total_dos: Optional[np.ndarray],
    grouped: Dict[PdosKey, np.ndarray],
    sp_counts: Dict[str, int],
    emin: float,
    emax: float,
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    mask = (E >= emin) & (E <= emax)
    Ewin = E[mask]

    max_x = 0.0
    if total_dos is not None:
        x = total_dos[mask]
        ax.fill_betweenx(Ewin, 0.0, x, color="red", alpha=0.18, linewidth=0.0, label="Total DOS (fill)")
        if x.size:
            max_x = max(max_x, float(np.nanmax(x)))

    keys = list(grouped.keys())
    color_by_key = build_color_map_for_keys(keys, sp_counts)
    for key in sorted(keys, key=lambda k: (k.species, k.inst, k.orbital)):
        y = grouped[key][mask]
        sp_lab = format_species_label(key.species, key.inst, sp_counts)
        ax.plot(y, Ewin, color=color_by_key[key], linewidth=1.3, alpha=0.95, label=f"{sp_lab} {key.orbital}")
        if y.size:
            max_x = max(max_x, float(np.nanmax(y)))

    if max_x <= 0.0 or (not np.isfinite(max_x)):
        max_x = 1.0

    ax.set_ylim(emin, emax)
    ax.set_xlim(0.0, 1.05 * max_x)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.set_xlabel("DOS / PDOS (arb.)")
    ax.set_ylabel("Energy (eV)")
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_components_for_each_file(
    outdir: Path,
    pdos_files: List[PdosFile],
    inst_map: Dict[Tuple[str, int], int],
    sp_counts: Dict[str, int],
    ef: Optional[float],
    align_fermi: bool,
    emin: float,
    emax: float,
    include_ldos: bool,
) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    saved = 0
    tab10 = plt.get_cmap("tab10")

    for pf in pdos_files:
        arr = _read_numeric_table(pf.path)
        if arr is None:
            continue

        E = arr[:, 0].copy()
        if align_fermi and (ef is not None):
            E = E - ef

        mask = (E >= emin) & (E <= emax)
        if not np.any(mask):
            continue

        inst = inst_map.get((pf.species, pf.atm), 1)
        sp_lab = format_species_label(pf.species, inst, sp_counts)

        fig, ax = plt.subplots(figsize=(6.6, 6.6))

        if include_ldos and arr.shape[1] >= 2:
            ax.plot(arr[mask, 1], E[mask], linewidth=1.0, color="0.35", alpha=0.9, label="ldos")

        ncol = arr.shape[1]
        if ncol >= 3:
            pdos_sum = np.nansum(arr[mask, 2:], axis=1)
            ax.fill_betweenx(E[mask], 0.0, pdos_sum, color="0.85", alpha=0.9, linewidth=0.0, label="pdos_sum (fill)")

            ncomp = ncol - 2
            labels = component_labels(ncomp=ncomp, orb_letter=pf.orb_letter)

            for j, col in enumerate(range(2, ncol)):
                ax.plot(
                    arr[mask, col],
                    E[mask],
                    color=tab10(j % 10),
                    linestyle="-",      # solid only
                    linewidth=1.35,
                    alpha=0.95,
                    label=labels[j] if j < len(labels) else f"pdos{col-1}",
                )

        ax.set_title(f"{sp_lab} wfc#{pf.wfc} ({pf.orb_raw})")
        ax.set_ylim(emin, emax)
        ax.set_xlim(left=0.0)
        ax.grid(True, linestyle=":", alpha=0.6)
        ax.set_xlabel("PDOS (arb.)")
        ax.set_ylabel("Energy (eV)")
        ax.legend(fontsize=9, frameon=True)

        outpath = outdir / f"{pf.path.name}_components.png"
        plt.tight_layout()
        plt.savefig(outpath, dpi=300)
        plt.close(fig)
        saved += 1

    return saved


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emin", type=float, default=None)
    ap.add_argument("--emax", type=float, default=None)

    ap.add_argument("--no-fermi-search", action="store_true")
    ap.add_argument("--fermi-from", default=None)
    ap.add_argument("--set-fermi", type=float, default=None)
    ap.add_argument("--no-align-fermi", action="store_true")

    ap.add_argument("--include-ldos", action="store_true")
    ap.add_argument("--overview-out", default="pdos_components/pdos_overview.png")
    ap.add_argument("--components-outdir", default="pdos_components/components")
    args = ap.parse_args()

    emin, emax = compute_ewin(args.emin, args.emax, default_emin=-3.0, default_emax=3.0)

    ef, ef_src = resolve_ef(args.no_fermi_search, args.fermi_from, args.set_fermi)
    align_fermi = (not args.no_align_fermi)

    prefix = find_prefix_from_pdos_tot(".")
    if prefix is None:
        print("No '*.pdos_tot' found in current directory. Cannot infer prefix.")
        return 2

    total_dos_path = find_total_dos_file(".")
    pdos_files = discover_all_projected_pdos(prefix, ".")
    if not pdos_files:
        print(f"No projected PDOS files found with prefix='{prefix}'.")
        return 2

    E, total_dos, grouped, sp_counts = load_grouped_pdos_on_grid(
        total_dos_path=total_dos_path,
        pdos_files=pdos_files,
        ef=ef,
        align_fermi=align_fermi,
    )
    if E is None:
        print("Failed to load any PDOS numeric data.")
        return 2

    plot_overview(Path(args.overview_out), E, total_dos, grouped, sp_counts, emin, emax)

    inst_map = renumber_species_instances(pdos_files)
    saved = plot_components_for_each_file(
        outdir=Path(args.components_outdir),
        pdos_files=pdos_files,
        inst_map=inst_map,
        sp_counts=sp_counts,
        ef=ef,
        align_fermi=align_fermi,
        emin=emin,
        emax=emax,
        include_ldos=args.include_ldos,
    )

    print(f"Ef source: {ef_src}, align_fermi={align_fermi}, Ef={ef}")
    print(f"Overview: {args.overview_out}")
    print(f"Components saved: {saved} -> {args.components_outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
