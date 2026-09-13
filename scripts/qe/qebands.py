#!/usr/bin/env python3
"""
Plot QE/plotband-style .gnu band file (segments) with:
- Default y-range = EF ± 2 eV when --yrange is not explicitly provided.
- Default x-range = full band-data extent [xmin, xmax] when --xrange is not set.
- X-axis uses HS labels as tick labels (numbers hidden).
- Automatic Fermi energy discovery from QE outputs (pw.x, dos.x, projwfc.x),
  preferring outputs whose filename shares the base prefix with the .gnu.
- Optional energy shift by EF (E -> E - EF) with informative message indicating source file.
- Automatic high-symmetry x-positions discovery from bands.x logs (vertical lines),
  and automatic HS labels from bands.in K_POINTS '!' comments or labelinfo.dat.
- Optional PDOS/DOS panel on the right when *.pdos_tot / *.pdos_atm#* files exist.

Usage
-----
python qebands.py band.dat.gnu \
  --yrange -12 12
# (auto-hs ON by default; bands.in / labelinfo auto-discovered)
"""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt


# ----------------------------
# I/O utilities
# ----------------------------
def read_text(path: Path) -> str:
    """Read file text with UTF-8 fallback to system encoding."""
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text()  # fallback


def split_segments(text: str) -> List[str]:
    """
    Split the file into segments by blank lines.
    A 'segment' corresponds to one continuous band chunk.
    """
    t = text.strip().replace("\r\n", "\n").replace("\r", "\n")
    return re.split(r"\n\s*\n+", t)


def clean_numeric_lines(segment: str) -> List[str]:
    """
    Keep only lines that start with numbers; drop comments/headers.
    Accepts leading spaces and +/-/decimal/exp formats.
    """
    out = []
    for line in segment.splitlines():
        s = line.strip()
        if not s:
            continue
        if re.match(r"^[\+\-]?\d+(\.\d+)?([eE][\+\-]?\d+)?\s+", s):
            out.append(s)
    return out


def load_segment_as_array(segment: str, xcol: int, ycol: int) -> np.ndarray:
    """
    Convert a cleaned numeric segment into an Nx2 array for chosen columns.
    xcol/ycol are 1-based indices.
    """
    data = np.genfromtxt(segment.splitlines())
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < max(xcol, ycol):
        raise ValueError("Segment has fewer columns than requested x/y columns.")
    return data[:, [xcol - 1, ycol - 1]]


def parse_gnu_file(path: Path, xcol: int, ycol: int) -> List[np.ndarray]:
    """Parse the .gnu file into a list of (x,y) arrays, one per segment."""
    text = read_text(path)
    segments = split_segments(text)
    arrays: List[np.ndarray] = []
    for seg in segments:
        numeric_lines = clean_numeric_lines(seg)
        if not numeric_lines:
            continue
        try:
            arr = load_segment_as_array("\n".join(numeric_lines), xcol, ycol)
            arr = arr[~np.isnan(arr).any(axis=1)]
            if arr.size:
                arrays.append(arr)
        except Exception:
            continue
    return arrays


# -----------------------------------------
# Discovery helpers
# -----------------------------------------
def _grep_first_float(pattern: re.Pattern, text: str) -> Optional[float]:
    """Return the first captured float from text or None."""
    m = pattern.search(text)
    if not m:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


class OutputCase:
    """Base interface for extracting values from a text output."""
    name: str = "base"
    def match(self, path: Path) -> bool: return False
    def extract(self, text: str) -> Optional[float]: return None


class CasePwLikeEF(OutputCase):
    """
    QE pw.x/dos.x/projwfc.x style:
    looks for: 'the Fermi energy is   <value> eV' (case-insensitive).
    """
    name = "pw/dos/projwfc:EF"
    PAT = re.compile(r"the\s+Fermi\s+energy\s+is\s+([\-+]?\d+(?:\.\d+)?)\s*eV", re.IGNORECASE)
    def match(self, path: Path) -> bool:
        s = path.name.lower()
        return any(k in s for k in ("scf.out", "nscf.out", "dos.out", "projwfc", ".out"))
    def extract(self, text: str) -> Optional[float]:
        return _grep_first_float(self.PAT, text)


CASES_EF: List[OutputCase] = [CasePwLikeEF()]


def _derive_base_prefix_from_gnu(gnu_file: Path) -> str:
    """
    Derive a base prefix from the .gnu filename to prefer matching .out/.in files.

    Examples
    --------
    mono_TiS2_4bands.dat.gnu -> mono_TiS2
    bulk_TiS2_bands.dat.gnu  -> bulk_TiS2
    sample.gnu                -> sample
    """
    name = gnu_file.name
    if name.endswith(".dat.gnu"):
        stem = name[:-len(".dat.gnu")]
    elif name.endswith(".gnu"):
        stem = name[:-len(".gnu")]
    else:
        stem = gnu_file.stem
    # Remove trailing _<num>bands, _bands, _<num>band
    stem = re.sub(r"(_\d*bands?)$", "", stem, flags=re.IGNORECASE)
    return stem


def discover_fermi_energy(search_dir: Path,
                          basename_hint: Optional[str] = None,
                          prefer_order: Sequence[str] = ("nscf", "scf", "dos", "projwfc", "out")
                          ) -> Tuple[Optional[float], Optional[Path]]:
    """
    Discover Fermi energy (in eV) by scanning outputs in `search_dir` with priority.

    Strategy
    --------
    1) If basename_hint is provided, prefer files whose name startswith that hint.
    2) Within each tier, apply prefer_order: nscf -> scf -> dos -> projwfc -> out.
    Returns:
      (EF, source_path) where EF is float or None, and source_path is the file used.
    """
    all_outs = sorted(search_dir.glob("*.out"))
    if not all_outs:
        return None, None

    def order_key(p: Path) -> Tuple[int, int]:
        ln = p.name.lower()
        tier = 1
        if basename_hint:
            if ln.startswith(basename_hint.lower()):
                tier = 0
        pr = 999
        for i, key in enumerate(prefer_order):
            if key in ln:
                pr = i
                break
        return (tier, pr)

    candidates = sorted(all_outs, key=lambda p: (*order_key(p), p.name.lower()))

    for p in candidates:
        txt = read_text(p)
        for case in CASES_EF:
            if case.match(p):
                ef = case.extract(txt)
                if ef is not None:
                    return ef, p
    return None, None


# ---- auto bands.in discovery ----
def discover_bands_in(search_dir: Path, basename_hint: Optional[str]) -> Optional[Path]:
    """
    Discover a plausible bands input (*.in) to read HS labels from K_POINTS.

    Priority:
      1) band.in, bands.in
      2) <basename_hint>*band*.in, <basename_hint>*.in (if hint provided)
      3) any *band*.in / *bands*.in in the current directory
      4) prefer files containing '1pw' when several candidates exist
      5) among candidates, prefer one that actually yields non-empty HS labels
    Returns:
      Path or None
    """
    candidates: List[Path] = []

    # 1) fixed names
    for fixed in ("band.in", "bands.in"):
        p = search_dir / fixed
        if p.exists():
            candidates.append(p)

    # 2) basename-hint-based candidates
    if basename_hint:
        candidates.extend(sorted(search_dir.glob(f"{basename_hint}*band*.in")))
        candidates.extend(sorted(search_dir.glob(f"{basename_hint}*.in")))

    # 3) generic candidates in current directory
    candidates.extend(sorted(search_dir.glob("*band*.in")))
    candidates.extend(sorted(search_dir.glob("*bands*.in")))

    # Deduplicate preserving order
    dedup: List[Path] = []
    seen = set()
    for p in candidates:
        if p.exists() and p not in seen:
            dedup.append(p)
            seen.add(p)

    if not dedup:
        return None

    def score(path: Path) -> Tuple[int, int, str]:
        name = path.name.lower()
        # smaller is better
        fixed_rank = 0 if name in ("band.in", "bands.in") else 1
        pw_rank = 0 if "1pw" in name else 1
        return (fixed_rank, pw_rank, name)

    # 4-5) try likely candidates first, but only accept files that actually parse labels
    for p in sorted(dedup, key=score):
        try:
            lbls = parse_bands_in_labels(p)
        except Exception:
            lbls = []
        if lbls:
            return p

    # If nothing produced labels, still return the best candidate for debugging/fallback
    return sorted(dedup, key=score)[0]


# ---- auto labelinfo discovery ----
def discover_labelinfo(search_dir: Path, basename_hint: Optional[str]) -> Optional[Path]:
    """Discover a plausible labelinfo file for HS labels."""
    candidates: List[Path] = []
    if basename_hint:
        candidates.extend(sorted(search_dir.glob(f"{basename_hint}*.labelinfo.dat")))
        candidates.extend(sorted(search_dir.glob(f"{basename_hint}*_band.labelinfo.dat")))
    candidates.extend(sorted(search_dir.glob("*_band.labelinfo.dat")))
    candidates.extend(sorted(search_dir.glob("*.labelinfo.dat")))

    seen = set()
    for path in candidates:
        if path.exists() and path not in seen:
            seen.add(path)
            return path
    return None


# ---- HS x from bands.x logs ----
class CaseBandsLogHS(OutputCase):
    """
    QE bands.x log style:
      'high-symmetry point: ... x coordinate   <value>'
    """
    name = "bands.x:HS"
    PAT = re.compile(r"high-symmetry\s+point:\s+.*?x\s+coordinate\s+([\-+]?\d+(?:\.\d+)?)",
                     re.IGNORECASE)
    def match(self, path: Path) -> bool:
        s = path.name.lower()
        return ("band" in s or "bands" in s) and s.endswith(".out")
    def extract_all(self, text: str) -> List[float]:
        return [float(v) for v in self.PAT.findall(text)]


CASES_HS: List[CaseBandsLogHS] = [CaseBandsLogHS()]


def discover_highsym_xcoords(search_dir: Path, bands_out: Optional[Path] = None) -> List[float]:
    """
    Discover high-symmetry x positions from bands.x logs.
    """
    paths: List[Path] = []
    if bands_out is not None:
        paths.append(bands_out)
    paths.extend(sorted(search_dir.glob("*band*.out")))
    paths.extend(sorted(search_dir.glob("*bands*.out")))
    dedup = []; seen = set()
    for p in paths:
        if p not in seen and p.exists():
            dedup.append(p); seen.add(p)
    xvals: List[float] = []
    for p in dedup:
        txt = read_text(p)
        for case in CASES_HS:
            if case.match(p):
                xs = case.extract_all(txt)
                if xs:
                    xvals.extend(xs)
    out: List[float] = []; seenf = set()
    for v in xvals:
        if v not in seenf:
            out.append(v); seenf.add(v)
    return out


# ----------------------------
# Parse HS labels from bands.in (K_POINTS)
# ----------------------------
def parse_bands_in_labels(bands_in: Path) -> List[str]:
    """
    Parse K_POINTS block from a QE bands input file and extract labels after '!'.

    Expected minimal pattern:
    K_POINTS <...>
      N
      kx ky kz w   ! Label
      ...
    """
    labels: List[str] = []
    txt = read_text(bands_in)

    m = re.search(r"^\s*K_POINTS\b.*$", txt, flags=re.IGNORECASE | re.MULTILINE)
    if not m:
        return labels

    start = m.end()
    lines = txt[start:].splitlines()

    # skip blanks/comments
    idx = 0
    while idx < len(lines) and (not lines[idx].strip() or lines[idx].lstrip().startswith(("!", "#"))):
        idx += 1
    if idx >= len(lines):
        return labels

    # read count if present
    n_nodes = None
    tok = lines[idx].strip().split()
    if tok and tok[0].isdigit():
        try:
            n_nodes = int(tok[0])
            idx += 1
        except Exception:
            n_nodes = None

    node_lines: List[str] = []
    if n_nodes is not None:
        for _ in range(n_nodes):
            if idx >= len(lines):
                break
            node_lines.append(lines[idx]); idx += 1
    else:
        # read until next card
        while idx < len(lines):
            line = lines[idx]
            if re.match(r"^\s*[A-Z_]+\b", line):
                break
            if line.strip():
                node_lines.append(line)
            idx += 1

    for nl in node_lines:
        if "!" in nl:
            lab = nl.split("!", 1)[1].strip()
            lab = re.sub(r"\s+", " ", lab)
            labels.append(lab if lab else "")
        else:
            labels.append("")
    return labels


# ----------------------------
# Parse HS labels from labelinfo.dat
# ----------------------------
@dataclass
class LabelInfo:
    label: str
    idx_1based: int
    kdist: float


def read_labelinfo_file(labelinfo_path: Path) -> List[LabelInfo]:
    """
    Expected format:
      LABEL   INDEX   KDIST   KX   KY   KZ ...
    """
    if not labelinfo_path.is_file():
        return []

    out: List[LabelInfo] = []
    with open(labelinfo_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                out.append(LabelInfo(label=parts[0], idx_1based=int(parts[1]), kdist=float(parts[2])))
            except Exception:
                continue
    return out


def labels_from_labelinfo(labelinfo_path: Path) -> List[str]:
    """Return labels in file order from labelinfo.dat."""
    return [item.label for item in read_labelinfo_file(labelinfo_path)]


# ----------------------------
# PDOS parsing / discovery
# ----------------------------
_PDOS_ATM_RE = re.compile(
    r"^(?P<prefix>.+?)\.pdos_atm\#(?P<atm>\d+)\((?P<species>[^)]+)\)_wfc\#(?P<wfc>\d+)\((?P<orb>[^)]+)\)$"
)


@dataclass(frozen=True)
class PdosKey:
    species: str
    species_idx: int
    orbital: str


@dataclass
class PdosData:
    energy: np.ndarray
    total: Optional[np.ndarray]
    projected: Dict[PdosKey, np.ndarray]


def _read_numeric_table(path: Path) -> Optional[np.ndarray]:
    rows: List[List[float]] = []
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#"):
                    continue
                parts = s.split()
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    continue
    except Exception:
        return None

    if not rows:
        return None
    arr = np.array(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        return None
    return arr


def find_total_pdos_file(search_dir: Path, basename_hint: Optional[str]) -> Optional[Path]:
    """Find *.pdos_tot, preferring basename-matched files when available."""
    cands = sorted(search_dir.glob("*.pdos_tot"))
    if not cands:
        return None
    if basename_hint:
        pref = [p for p in cands if p.name.startswith(basename_hint)]
        if pref:
            return pref[0]
    return cands[0]


def discover_projected_pdos_files(search_dir: Path, prefix: str) -> List[Tuple[Path, int, str, int, str]]:
    """
    Return list of tuples:
      (path, atm_index_original, species, wfc_index, orbital_letter)
    """
    out: List[Tuple[Path, int, str, int, str]] = []
    for p in sorted(search_dir.glob(f"{prefix}.pdos_atm#*")):
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
    by_species: Dict[str, List[int]] = {}
    for _, atm, species, _, _ in entries:
        by_species.setdefault(species, []).append(atm)

    mapping: Dict[Tuple[str, int], int] = {}
    for species, atms in by_species.items():
        for idx, atm in enumerate(sorted(set(atms)), start=1):
            mapping[(species, atm)] = idx
    return mapping


def species_to_cmap_name(species_list: List[str]) -> Dict[str, str]:
    palette = ["spring", "summer", "winter", "autumn", "cool", "copper"]
    out: Dict[str, str] = {}
    for i, species in enumerate(species_list):
        out[species] = palette[i % len(palette)]
    return out


def load_pdos_grouped(search_dir: Path,
                      basename_hint: Optional[str],
                      ef: Optional[float],
                      shift_by_ef: bool) -> Optional[PdosData]:
    total_pdos_path = find_total_pdos_file(search_dir, basename_hint)
    proj_entries: List[Tuple[Path, int, str, int, str]] = []
    prefix = None

    if total_pdos_path is not None:
        prefix = total_pdos_path.name[:-len(".pdos_tot")]
        proj_entries = discover_projected_pdos_files(search_dir, prefix)

    if total_pdos_path is None and not proj_entries:
        return None

    energy_grid: Optional[np.ndarray] = None
    total_dos: Optional[np.ndarray] = None

    if total_pdos_path is not None:
        arr = _read_numeric_table(total_pdos_path)
        if arr is not None:
            energy_grid = arr[:, 0].copy()
            total_dos = arr[:, 1].copy()

    if energy_grid is None and proj_entries:
        arr0 = _read_numeric_table(proj_entries[0][0])
        if arr0 is not None:
            energy_grid = arr0[:, 0].copy()

    if energy_grid is None:
        return None

    if shift_by_ef and ef is not None:
        energy_grid = energy_grid - ef

    projected: Dict[PdosKey, np.ndarray] = {}
    inst_map = renumber_species_instances(proj_entries)

    for path, atm, species, _wfc, orb_letter in proj_entries:
        arr = _read_numeric_table(path)
        if arr is None:
            continue
        e = arr[:, 0].copy()
        y = np.nansum(arr[:, 1:], axis=1)
        if shift_by_ef and ef is not None:
            e = e - ef

        order = np.argsort(e)
        e_sorted = e[order]
        y_sorted = y[order]
        e_unique, idx_unique = np.unique(e_sorted, return_index=True)
        y_unique = y_sorted[idx_unique]
        y_grid = np.interp(energy_grid, e_unique, y_unique, left=0.0, right=0.0)

        key = PdosKey(species=species,
                      species_idx=inst_map.get((species, atm), 1),
                      orbital=orb_letter)
        if key in projected:
            projected[key] = projected[key] + y_grid
        else:
            projected[key] = y_grid

    if total_dos is None and not projected:
        return None

    return PdosData(energy=energy_grid, total=total_dos, projected=projected)


# ----------------------------
# Plotting (HS labels as x-tick labels)
# ----------------------------
def plot_segments(
    arrays: List[np.ndarray],
    yrange: Tuple[float, float] = None,
    xrange: Tuple[float, float] = None,
    color: str = "red",
    lw: float = 2.5,
    xlabel: str = "",
    ylabel: str = "Energy (eV)",
    title: str = "Band structure (segments as continuous curves)",
    show_grid: bool = True,
    hs_x: Optional[List[float]] = None,
    hs_labels: Optional[List[str]] = None,
    hide_x_axis_numbers_when_no_hs: bool = True,
    pdos_data: Optional[PdosData] = None,
) -> plt.Figure:
    """
    Make a matplotlib figure plotting each segment as a separate continuous line.
    - If hs_x is provided, use hs_x as xticks and hs_labels as tick labels.
    - No in-plot HS text; only x-axis tick labels.
    - If pdos_data is provided, add a DOS/PDOS panel on the right.
    """
    if pdos_data is not None:
        fig = plt.figure(figsize=(9, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        ax_dos = fig.add_subplot(gs[0, 1], sharey=ax)
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax_dos = None

    for arr in arrays:
        x, y = arr[:, 0], arr[:, 1]
        ax.plot(x, y, "-", linewidth=lw, color=color)

    # x-range
    if xrange is not None:
        ax.set_xlim(*xrange)

    # y-range
    if yrange is not None:
        ax.set_ylim(*yrange)

    ax.tick_params(axis='y', labelsize=14)
    ax.yaxis.set_major_locator(plt.MaxNLocator(6, symmetric=True, integer=True))

    # HS ticks and labels
    if hs_x and len(hs_x) >= 1:
        if not hs_labels or len(hs_labels) != len(hs_x):
            hs_labels = [f"HS{i+1}" for i in range(len(hs_x))]
        ax.set_xticks(hs_x)
        ax.set_xticklabels(hs_labels, fontsize=20)
        ax.tick_params(axis="x", which="both", length=4)
        ax.set_xlabel("")
    else:
        if hide_x_axis_numbers_when_no_hs:
            ax.set_xticks([])
            ax.tick_params(axis="x", which="both", length=0)
        ax.set_xlabel("")

    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_facecolor("none")
    if show_grid:
        ax.grid(alpha=0.3)

    if ax_dos is not None and pdos_data is not None:
        max_x = 0.0
        if pdos_data.total is not None:
            ax_dos.plot(pdos_data.total, pdos_data.energy, color="red", linewidth=1.6, alpha=0.95, label="Total DOS")
            try:
                max_x = max(max_x, float(np.nanmax(pdos_data.total)))
            except Exception:
                pass

        species_order = sorted({key.species for key in pdos_data.projected.keys()})
        cmap_by_species = species_to_cmap_name(species_order)

        by_inst: Dict[Tuple[str, int], List[PdosKey]] = {}
        for key in pdos_data.projected.keys():
            by_inst.setdefault((key.species, key.species_idx), []).append(key)

        for species, species_idx in sorted(by_inst.keys(), key=lambda item: (item[0], item[1])):
            keys = sorted(by_inst[(species, species_idx)], key=lambda k: k.orbital)
            cmap = plt.get_cmap(cmap_by_species.get(species, "viridis"))
            sample_points = np.linspace(0.15, 0.95, max(1, len(keys)))
            for j, key in enumerate(keys):
                y = pdos_data.projected[key]
                ax_dos.plot(y, pdos_data.energy,
                            color=cmap(sample_points[j]),
                            linewidth=1.2,
                            alpha=0.95,
                            label=f"{key.species}#{key.species_idx} {key.orbital}")
                try:
                    max_x = max(max_x, float(np.nanmax(y)))
                except Exception:
                    pass

        if max_x <= 0.0 or not np.isfinite(max_x):
            max_x = 1.0
        ax_dos.set_xlim(0.0, 1.05 * max_x)
        ax_dos.grid(alpha=0.3)
        ax_dos.set_xlabel("DOS / PDOS", fontsize=14)
        plt.setp(ax_dos.get_yticklabels(), visible=False)
        ax_dos.tick_params(axis="x", labelsize=10)
        if pdos_data.total is not None or pdos_data.projected:
            ax_dos.legend(loc="upper right", fontsize=8, frameon=True, borderpad=0.3, handlelength=1.5)

    if pdos_data is None:
        fig.tight_layout()
    else:
        fig.subplots_adjust(left=0.10, right=0.97, bottom=0.12, top=0.98, wspace=0.05)
    return fig


# ----------------------------
# Output name derivation
# ----------------------------
def _derive_outputs(gnu_file: Path) -> tuple[Path, Path]:
    """
    Derive PNG/PDF output paths from the input .gnu-like filename.
    """
    name = gnu_file.name
    parent = gnu_file.parent
    if name.endswith(".dat.gnu"):
        stem = name[:-len(".dat.gnu")]
    elif name.endswith(".gnu"):
        stem = name[:-len(".gnu")]
    else:
        stem = gnu_file.stem
    return parent / f"{stem}.png", parent / f"{stem}.pdf"


# ----------------------------
# CLI
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Plot .gnu band file without inter-segment lines; EF-centered y-range, HS vlines & HS tick labels."
    )
    p.add_argument("gnu_file", nargs="?", type=Path, default=None,
                   help="Input .gnu file (default: auto-detect band.dat.gnu, then *.gnu in current directory)")
    p.add_argument("--xcol", type=int, default=1, help="1-based column index for x (default: 1)")
    p.add_argument("--ycol", type=int, default=2, help="1-based column index for y (default: 2)")

    # ranges
    p.add_argument("--yrange", nargs=2, type=float, metavar=("YMIN", "YMAX"),
                   help="Set y-range; if omitted, defaults to EF±3 eV when EF is found.")
    p.add_argument("--xrange", nargs=2, type=float, metavar=("XMIN", "XMAX"),
                   help="Set x-range explicitly (e.g., --xrange 0.0 3.5).")

    # style
    p.add_argument("--color", default="red", help="Line color (default: red)")
    p.add_argument("--lw", type=float, default=2.5, help="Line width (default: 2.5)")

    # outputs (derived from input if omitted)
    p.add_argument("--png", type=Path, default=None, help="PNG output path (default: derive from input)")
    p.add_argument("--pdf", type=Path, default=None, help="PDF output path (default: derive from input)")

    # EF discovery/shift
    p.add_argument("--ef-auto", default=True, action="store_true",
                   help="Prefer to discover Fermi energy from QE outputs in the same directory.")
    p.add_argument("--ef-source", type=Path,
                   help="Explicit output file to read EF from (overrides auto discovery).")
    p.add_argument("--shift-by-ef", dest="shift_by_ef", default=True, action="store_true",
                   help="Shift energy by EF (E -> E - EF) after discovery.")
    p.add_argument("--no-shift-by-ef", dest="shift_by_ef", action="store_false",
                   help="Do not shift energy by EF even if EF is found.")

    # HS discovery (default ON)
    p.add_argument("--auto-hs", dest="auto_hs", action="store_true", default=True,
                   help="Auto-discover high-symmetry x positions from bands.x logs and draw vlines (default ON).")
    p.add_argument("--no-auto-hs", dest="auto_hs", action="store_false",
                   help="Disable automatic high-symmetry discovery.")
    p.add_argument("--bands-out", type=Path,
                   help="Explicit bands.x output file to parse for high-symmetry x positions.")
    p.add_argument("--bands-in", type=Path,
                   help="QE bands input file containing K_POINTS node labels (comments after '!').")
    p.add_argument("--labelinfo", type=Path,
                   help="Explicit labelinfo.dat file containing HS labels.")
    p.add_argument("--hs-labels", type=str,
                   help='Comma-separated labels for HS points, e.g., "Γ,M,K,Γ".')

    args = p.parse_args()

    # ---- .gnu input auto-discovery when omitted ----
    if args.gnu_file is None:
        preferred = Path.cwd() / "band.dat.gnu"
        if preferred.is_file():
            args.gnu_file = preferred
        else:
            gnu_candidates = sorted(Path.cwd().glob("*.gnu"))
            if not gnu_candidates:
                raise SystemExit("No .gnu file found in current directory.")
            args.gnu_file = gnu_candidates[0]
        print(f"[INFO] Auto-detected .gnu file: {args.gnu_file.name}")

    arrays = parse_gnu_file(args.gnu_file, xcol=args.xcol, ycol=args.ycol)
    if not arrays:
        raise SystemExit("No numeric segments found. Check your input file or column indices.")

    # ---- EF discovery ----
    ef_val: Optional[float] = None
    ef_src: Optional[Path] = None
    if args.ef_source:
        txt = read_text(args.ef_source)
        for case in CASES_EF:
            if case.match(args.ef_source):
                ef_val = case.extract(txt); ef_src = args.ef_source if ef_val is not None else None
                break
    else:
        base_prefix = _derive_base_prefix_from_gnu(args.gnu_file)
        ef_val, ef_src = discover_fermi_energy(args.gnu_file.parent, basename_hint=base_prefix)

    # ---- EF shift (optional) ----
    if ef_val is not None and args.shift_by_ef:
        shifted = []
        for arr in arrays:
            arr2 = arr.copy()
            arr2[:, 1] = arr2[:, 1] - ef_val
            shifted.append(arr2)
        arrays = shifted
        src_name = ef_src.name if ef_src is not None else "unknown"
        print(f"[INFO] Applied EF shift: E -> E - {ef_val:.6f} eV (found at {src_name})")
    elif ef_val is None and args.shift_by_ef:
        print("[WARN] --shift-by-ef requested but EF not found; skipping shift.")

    # ---- y-range defaulting to EF ± 2 eV if not specified ----
    chosen_yrange: Optional[Tuple[float, float]] = None
    if args.yrange is not None:
        chosen_yrange = tuple(args.yrange)
    else:
        if ef_val is not None:
            if args.shift_by_ef:
                chosen_yrange = (-3.0, 3.0)
            else:
                chosen_yrange = (ef_val - 3.0, ef_val + 3.0)

    # ---- HS discovery (x coords + labels) ----
    hs_x: Optional[List[float]] = None
    hs_labels: Optional[List[str]] = None
    bands_in_used: Optional[Path] = None
    labelinfo_used: Optional[Path] = None

    if args.auto_hs:
        hs_x = discover_highsym_xcoords(args.gnu_file.parent, args.bands_out)
        if not hs_x:
            print("[WARN] No high-symmetry x positions found from bands.x logs.")

        if args.hs_labels:
            parsed = [s.strip() for s in args.hs_labels.split(",")]
            if hs_x and len(parsed) != len(hs_x):
                print("[WARN] --hs-labels count does not match HS points; ignoring custom labels.")
            else:
                hs_labels = parsed

        if hs_labels is None:
            labelinfo_file = args.labelinfo
            if labelinfo_file is None:
                base_prefix = _derive_base_prefix_from_gnu(args.gnu_file)
                labelinfo_file = discover_labelinfo(args.gnu_file.parent, base_prefix)
            if labelinfo_file and labelinfo_file.exists():
                lbls = labels_from_labelinfo(labelinfo_file)
                labelinfo_used = labelinfo_file
                if hs_x and len(lbls) != len(hs_x):
                    print("[WARN] labelinfo labels count does not match HS points; falling back to bands.in / auto labels.")
                else:
                    hs_labels = lbls if lbls else None

        if hs_labels is None:
            bands_in_file = args.bands_in
            if bands_in_file is None:
                base_prefix = _derive_base_prefix_from_gnu(args.gnu_file)
                bands_in_file = discover_bands_in(args.gnu_file.parent, base_prefix)
            if bands_in_file and bands_in_file.exists():
                lbls = parse_bands_in_labels(bands_in_file)
                bands_in_used = bands_in_file
                if hs_x and len(lbls) != len(hs_x):
                    print("[WARN] bands.in labels count does not match HS points; falling back to auto labels.")
                else:
                    hs_labels = lbls if lbls else None

    # ---- default x-range = full band-data extent if not provided ----
    chosen_xrange: Optional[Tuple[float, float]] = None
    if args.xrange is not None:
        chosen_xrange = tuple(args.xrange)
    else:
        xmins = [float(np.nanmin(arr[:, 0])) for arr in arrays if arr.size]
        xmaxs = [float(np.nanmax(arr[:, 0])) for arr in arrays if arr.size]
        if xmins and xmaxs:
            chosen_xrange = (min(xmins), max(xmaxs))

    # ---- PDOS auto-discovery ----
    base_prefix = _derive_base_prefix_from_gnu(args.gnu_file)
    pdos_data = load_pdos_grouped(args.gnu_file.parent, base_prefix, ef=ef_val, shift_by_ef=bool(args.shift_by_ef))
    total_pdos_file = find_total_pdos_file(args.gnu_file.parent, base_prefix)
    projected_pdos_files: List[Path] = []
    if total_pdos_file is not None:
        projected_pdos_files = [item[0] for item in discover_projected_pdos_files(
            args.gnu_file.parent, total_pdos_file.name[:-len('.pdos_tot')]
        )]

    # ---- derive output paths if omitted ----
    png_path, pdf_path = (args.png, args.pdf)
    if png_path is None or pdf_path is None:
        d_png, d_pdf = _derive_outputs(args.gnu_file)
        if png_path is None: png_path = d_png
        if pdf_path is None: pdf_path = d_pdf

    # ---- plot ----
    fig = plot_segments(
        arrays,
        yrange=chosen_yrange,
        xrange=chosen_xrange,
        color=args.color,
        lw=args.lw,
        hs_x=hs_x,
        hs_labels=hs_labels,
        hide_x_axis_numbers_when_no_hs=True,
        pdos_data=pdos_data,
    )

    # optional info prints
    if args.auto_hs and hs_x:
        print(f"[INFO] High-symmetry x positions: {', '.join(f'{v:.4f}' for v in hs_x)}")
        if hs_labels:
            print(f"[INFO] High-symmetry labels: {', '.join(hs_labels)}")
        if labelinfo_used:
            print(f"[INFO] HS labels source: {labelinfo_used.name}")
        elif bands_in_used:
            print(f"[INFO] HS labels source: {bands_in_used.name}")

    if total_pdos_file is not None or projected_pdos_files:
        tot_name = total_pdos_file.name if total_pdos_file is not None else 'None'
        print(f"[INFO] PDOS autodetect: total={tot_name}, projected={len(projected_pdos_files)} file(s)")
        if projected_pdos_files:
            preview = ', '.join(path.name for path in projected_pdos_files[:4])
            extra = '' if len(projected_pdos_files) <= 4 else ', ...'
            print(f"[INFO] PDOS projected files: {preview}{extra}")
    else:
        print("[INFO] PDOS autodetect: no *.pdos_tot / *.pdos_atm#* files found")

    fig.patch.set_alpha(0.0)
    fig.savefig(png_path, dpi=300, transparent=True)
    fig.savefig(pdf_path)
    print(f"Saved: {png_path} and {pdf_path}")


if __name__ == "__main__":
    main()


