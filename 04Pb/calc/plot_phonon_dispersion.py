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
        "--freq", default=None,
        help=(
            "Path to .freq.gp (preferred) or raw .freq file, or PREFIX. "
            "If omitted, all *.freq.gp files in the current directory are listed for selection; "
            "after 10 seconds with no response, the first item is selected automatically."
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


def _infer_prefix_from_freq_path(freq_path: Path) -> str:
    """Infer a QE prefix from common dispersion filenames."""
    name = freq_path.name
    suffixes = [
        "_phband.freq.gp.dat", "_phband.freq.dat",
        "_phband.freq.gp", "_phband.freq",
        "_dispersion.freq.gp", "_dispersion.freq",
        ".freq.gp", ".freq",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return freq_path.stem


def choose_freq_gp_interactively(timeout_sec: int = 10) -> Path:
    """
    If --freq is omitted, list *.freq.gp files in the current directory and
    let the user choose one. If no input arrives within timeout_sec, select
    the first item automatically.
    """
    import signal

    candidates = sorted(Path.cwd().glob("*.freq.gp"), key=lambda p: p.name.lower())
    if not candidates:
        raise FileNotFoundError(
            "--freq was not given and no *.freq.gp file was found in the current directory."
        )

    print("[INFO] Found *.freq.gp files:")
    for i, p in enumerate(candidates, start=1):
        print(f"  {i:2d}) {p.name}")

    default = candidates[0]
    prompt = (
        f"Select dispersion file [1-{len(candidates)}] "
        f"(default: 1 = {default.name}, timeout: {timeout_sec}s): "
    )

    class _InputTimeout(Exception):
        pass

    def _alarm_handler(signum, frame):
        raise _InputTimeout

    answer = ""
    old_handler = None
    try:
        # Linux/HPC use case: SIGALRM gives a real timeout around input().
        old_handler = signal.signal(signal.SIGALRM, _alarm_handler)
        signal.alarm(timeout_sec)
        answer = input(prompt).strip()
        signal.alarm(0)
    except (_InputTimeout, EOFError, KeyboardInterrupt):
        signal.alarm(0)
        print(f"\n[INFO] No selection -> using {default.name}")
        return default
    finally:
        if old_handler is not None:
            signal.signal(signal.SIGALRM, old_handler)

    if not answer:
        print(f"[INFO] Empty input -> using {default.name}")
        return default

    try:
        idx = int(answer)
        if 1 <= idx <= len(candidates):
            chosen = candidates[idx - 1]
            print(f"[INFO] Selected: {chosen.name}")
            return chosen
    except ValueError:
        pass

    print(f"[WARN] Invalid selection '{answer}' -> using {default.name}")
    return default


def _strip_inline_comment(line: str) -> Tuple[str, str]:
    """Split a QE-style line into data and an optional ! comment."""
    if "!" in line:
        data, comment = line.split("!", 1)
        return data.strip(), comment.strip()
    return line.strip(), ""


def _normalize_hs_label(label: str, index: int) -> str:
    """Normalize common Gamma spellings while preserving user labels."""
    label = label.strip()
    if not label:
        return f"P{index}"
    if re.fullmatch(r"(?i)(g|gamma|Γ)", label):
        return "G"
    return label


def parse_matdyn_band_qpath(path: Path) -> Optional[List[Tuple[float, float, float, int, str]]]:
    """
    Parse a QE matdyn/ph-style input containing q_in_band_form=.true. and a
    post-namelist q-point block:

        /
        N
        qx qy qz npts ! LABEL
    """
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None

    text = "\n".join(lines)
    if not re.search(r"q_in_band_form\s*=\s*\.true\.", text, flags=re.I):
        return None

    slash_idx = next((i for i, raw in enumerate(lines) if raw.strip() == "/"), None)
    if slash_idx is None:
        return None

    j = slash_idx + 1
    while j < len(lines):
        data, _ = _strip_inline_comment(lines[j])
        if data and not data.startswith("#"):
            break
        j += 1
    if j >= len(lines):
        return None

    try:
        n_hsp = int(data.split()[0])
    except (ValueError, IndexError):
        return None

    rows: List[Tuple[float, float, float, int, str]] = []
    j += 1
    while j < len(lines) and len(rows) < n_hsp:
        raw = lines[j].strip()
        j += 1
        if not raw or raw.startswith("#"):
            continue

        data, comment = _strip_inline_comment(raw)
        toks = data.split()
        if len(toks) < 4:
            continue
        try:
            qx, qy, qz = map(float, toks[:3])
            npts = int(float(toks[3]))
        except ValueError:
            return None

        label = comment.split()[0] if comment else (toks[4] if len(toks) >= 5 else "")
        label = _normalize_hs_label(label, len(rows) + 1)
        rows.append((qx, qy, qz, npts, label))

    return rows if len(rows) == n_hsp and n_hsp >= 2 else None


def discover_qpath_source() -> Optional[Tuple[Path, List[Tuple[float, float, float, int, str]]]]:
    """Search current directory for parseable files whose names contain ph or matdyn."""
    candidates = [
        p for p in Path.cwd().iterdir()
        if p.is_file() and ("ph" in p.name.lower() or "matdyn" in p.name.lower())
    ]
    candidates.sort(
        key=lambda p: (
            0 if p.suffix.lower() == ".in" else 1,
            0 if "matdyn" in p.name.lower() else 1,
            p.name.lower(),
        )
    )

    for p in candidates:
        rows = parse_matdyn_band_qpath(p)
        if rows:
            return p, rows
    return None


def write_qpath_in(
    rows: List[Tuple[float, float, float, int, str]],
    out_path: Path = Path("qpath.in"),
) -> Path:
    """Write qpath.in for bandplot / plot_phonon_dispersion.py."""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write('"""\n')
        f.write("qpath.in for bandplot & plot_phonon_dispersion.py\n")
        f.write('"""\n')
        f.write(f"{len(rows)}\n")
        for qx, qy, qz, npts, label in rows:
            f.write(f"{qx:.6f}  {qy:.6f}  {qz:.6f}  {npts:d}  {label}\n")
    return out_path


def ensure_qpath_in(qpath_arg: Optional[str]) -> Optional[Path]:
    """Resolve --qpath, or create ./qpath.in from local *ph*/*matdyn* input."""
    if qpath_arg:
        p = Path(qpath_arg)
        if p.exists() and p.is_file():
            return p
        print(f"[WARN] --qpath given but not found. Skipping HS labels: {qpath_arg}")
        return None

    default = Path("qpath.in")
    if default.exists() and default.is_file():
        return default

    found = discover_qpath_source()
    if found is None:
        print(
            "[WARN] qpath.in not found, and no parseable band-form q-path was found "
            "in local files containing 'ph' or 'matdyn'. Skipping HS labels."
        )
        return None

    src, rows = found
    out = write_qpath_in(rows, default)
    print(f"[OK] Generated {out} from {src.name}")
    return out


def resolve_inputs(
    freq_arg: Optional[str], dos_arg: Optional[str], qpath_arg: Optional[str]
) -> Tuple[Path, Optional[Path], Optional[Path], str]:
    """Resolve dispersion/DOS/qpath files from the current working directory."""
    prefix = ""

    if freq_arg is None:
        freq_path = choose_freq_gp_interactively()
        prefix = _infer_prefix_from_freq_path(freq_path)
    else:
        freq_candidate = Path(freq_arg)
        if freq_candidate.exists() and freq_candidate.is_file():
            freq_path = freq_candidate
            prefix = _infer_prefix_from_freq_path(freq_path)
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

    dos_path: Optional[Path] = None
    if dos_arg:
        p = Path(dos_arg)
        if p.exists() and p.is_file():
            dos_path = p
        else:
            print(f"[WARN] --dos given but not found. Skipping DOS: {dos_arg}")
    elif prefix:
        dos_guesses = [
            Path(f"{prefix}_phdos"),
            Path(f"{prefix}.phdos"),
            Path(f"{prefix}_phdos.dat"),
            Path(f"{prefix}.phdos.dat"),
        ]
        dos_path = _first_existing(dos_guesses)

    qpath_path = ensure_qpath_in(qpath_arg)
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
            raw_lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

        # qpath.in may begin with a small triple-quoted description block.
        lines: List[str] = []
        in_doc = False
        for ln in raw_lines:
            if ln.startswith('"""'):
                in_doc = not in_doc
                continue
            if not in_doc:
                lines.append(ln)

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
def _style_dispersion_axis(ax, q_path, labels, tick_idx, ylabel, full_title):
    """Apply the reference-paper style used for the phonon dispersion panel."""
    ax.set_xlim(float(np.min(q_path)), float(np.max(q_path)))
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_title(full_title, fontsize=24, pad=12)

    # Large, clean ticks; HS labels are intentionally prominent.
    ax.tick_params(axis="y", labelsize=16, width=1.2, length=5)
    ax.tick_params(axis="x", labelsize=24, width=1.2, length=0, pad=6)

    # Light horizontal dashed guides only, matching the reference figure.
    ax.grid(axis="y", linestyle="--", linewidth=0.8, color="0.82", alpha=0.9)
    ax.grid(axis="x", visible=False)

    # Darker vertical high-symmetry separators.
    if labels and tick_idx and len(labels) == len(tick_idx):
        tick_pos = [q_path[i] for i in tick_idx]
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(labels, fontsize=24)
        for x in tick_pos:
            ax.axvline(x, linestyle="-", linewidth=1.25, color="0.45", alpha=0.95, zorder=0)
    else:
        ax.set_xticks([])

    # Strong rectangular frame.
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("black")


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

    # Reference style: thin vivid-red phonon branches, no markers.
    for j in range(y_disp.shape[1]):
        ax.plot(q_path, y_disp[:, j], "-", lw=1.6, color="red", solid_capstyle="round")
    ax.axhline(0.0, lw=1.0, color="0.35", zorder=0)

    if (emin is not None) or (emax is not None):
        lo = emin if emin is not None else float(np.nanmin(y_disp))
        hi = emax if emax is not None else float(np.nanmax(y_disp))
        ax.set_ylim(lo, hi)

    _style_dispersion_axis(ax, q_path, labels, tick_idx, ylabel, full_title)
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
        ax.plot(q_path, y_disp[:, j], "-", lw=1.6, color="red", solid_capstyle="round")
    ax.axhline(0.0, lw=1.0, color="0.35", zorder=0)

    if (emin is not None) or (emax is not None):
        lo = emin if emin is not None else float(np.nanmin(y_disp))
        hi = emax if emax is not None else float(np.nanmax(y_disp))
        ax.set_ylim(lo, hi)

    _style_dispersion_axis(ax, q_path, labels, tick_idx, ylabel, full_title)

    # DOS panel: same red line language, restrained fill and matching frame.
    ax_dos.plot(dos_val, y_dos, "-", lw=1.6, color="red", solid_capstyle="round")
    ax_dos.fill_betweenx(y_dos, 0.0, dos_val, alpha=0.10, color="red")
    ax_dos.set_xlabel("DOS", fontsize=16)
    ax_dos.grid(axis="y", linestyle="--", linewidth=0.8, color="0.82", alpha=0.9)
    ax_dos.grid(axis="x", visible=False)
    ax_dos.tick_params(axis="x", labelsize=12, width=1.2, length=4)
    plt.setp(ax_dos.get_yticklabels(), visible=False)
    ax_dos.tick_params(axis="y", length=0)
    for spine in ax_dos.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("black")

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
