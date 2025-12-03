#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QE PDOS overlay plotter with automatic Fermi-level shift and spin handling.

Features
--------
1. Scan current directory for QE PDOS files with standard projwfc.x naming:
   - prefix.pdos_tot[_up/down]
   - prefix.pdos_atm#N(Element)_wfc#M(orbital)[_up/down]
   Any other "pdos" files (e.g. prefix.pdos.projwfc_up, *.in, *.out, *.py)
   are ignored.

2. Detect Fermi energy from QE output files (*.out):
   - Search for lines: "the Fermi energy is   X eV".
   - If found, optionally shift E -> E - EF based on a heuristic:
        if energy grid already looks symmetric around 0, do NOT shift.

3. Spin handling:
   - If no spin-resolved PDOS files (no *_up / *_down suffix):
       → plot normal DOS: all curves positive.
   - If spin-resolved PDOS files exist:
       → treat as spin-polarized:
          * spin up:  y(E)  > 0
          * spin down: y(E) < 0 (we plot -DOS for ↓)
       → y-axis label explicitly mentions spin-down plotted negative.

4. Single figure with:
   - Total DOS (if present)
   - All atom/orbital PDOS curves overlaid
   - Energy axis either raw E or (E - EF), with deterministic x-label.

5. Output filename:
   - Automatically uses the detected seedname from PDOS files:
       seedname_pdos_overlay.png

Limitations
-----------
- For projwfc.x PDOS files, we only use the first DOS column after E.
  If you want to sum over all m-resolved states, modify the lines:
      y = cols[:, 0]
  to something like:
      y = cols.sum(axis=1)
  depending on your preferred definition.

- The EF-shift decision is heuristic. Always visually check whether
  0 eV is indeed the Fermi level.
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Parsing utilities
# -----------------------------

def parse_pdos_file(filepath):
    """
    Read a QE PDOS file and return (E, DOS_columns, header_line).

    Parameters
    ----------
    filepath : str
        Path to the PDOS file.

    Returns
    -------
    E : np.ndarray, shape (N,)
        Energy grid in eV (as written in the file).
    cols : np.ndarray, shape (N, M)
        DOS columns (all columns except energy).
    header : str or None
        Last header line starting with '#', if found; otherwise None.
    """
    data = []
    header = None
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                header = line.strip()
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                row = list(map(float, parts))
            except ValueError:
                # Non-numeric lines are silently skipped
                continue
            data.append(row)

    if not data:
        raise ValueError(f"No numerical data in file {filepath}")

    arr = np.array(data, dtype=float)
    return arr[:, 0], arr[:, 1:], header


def extract_label_from_filename(filename):
    """
    Generate a short label for PDOS curves based on QE filename.

    Examples
    --------
    'TiS2_monolayer.pdos_atm#1(Ti)_wfc#3(d)'      -> 'Ti d'
    'TiS2_monolayer.pdos_atm#2(S)_wfc#2(p)_up'   -> 'S p'
    'TiS2_monolayer.pdos_tot'                    -> 'Total DOS'
    'TiS2_monolayer.pdos_tot_up'                 -> 'Total DOS'
    """
    # Total DOS
    if ".pdos_tot" in filename:
        return "Total DOS"

    # Atom symbol inside parentheses: (Ti), (S), ...
    atm = re.search(r"\((.*?)\)", filename)  # first '(...)'
    atom = atm.group(1) if atm else "atom?"

    # Orbital from 'wfc#n(l)'
    orb = re.search(r"wfc#\d+\((.*?)\)", filename)
    orbital = orb.group(1) if orb else "wfc?"

    return f"{atom} {orbital}"


# -----------------------------
# Fermi energy detection
# -----------------------------

def detect_fermi_energy_from_out():
    """
    Search QE output files (*.out) for a line containing 'the Fermi energy is'.

    Returns
    -------
    ef : float or None
        Detected Fermi energy in eV, or None if not found.

    Notes
    -----
    - Assumes QE's standard printout:
         'the Fermi energy is   XX.XXXX eV'
    - If multiple output files contain a Fermi energy, take the last
      occurrence found (usually the most recent run).
    """
    pattern = re.compile(r"the Fermi energy is\s+([0-9Ee\+\-\.]+)\s+eV")
    ef_values = []

    for fname in os.listdir("."):
        if not fname.endswith(".out"):
            continue
        try:
            with open(fname, "r") as f:
                for line in f:
                    m = pattern.search(line)
                    if m:
                        try:
                            ef = float(m.group(1))
                            ef_values.append(ef)
                        except ValueError:
                            continue
        except (IOError, OSError):
            continue

    if ef_values:
        return ef_values[-1]
    return None


def should_apply_fermi_shift(E, ef):
    """
    Decide heuristically whether we should apply a Fermi shift.

    Parameters
    ----------
    E : np.ndarray
        Raw energy grid from a PDOS file.
    ef : float
        Fermi energy in eV obtained from QE outputs.

    Returns
    -------
    apply_shift : bool
        True if we decide to use E_shifted = E - ef, False otherwise.

    Heuristic (explicitly a guess):
    - If the energy grid spans significant negative and positive energies
      and is roughly symmetric around 0, we assume PDOS is already EF-centered.
    - Otherwise, we decide to apply the shift.
    """
    Emin = float(E.min())
    Emax = float(E.max())

    if Emin < 0.0 < Emax:
        left = abs(Emin)
        right = abs(Emax)
        if left > 1e-6 and right > 1e-6:
            ratio = max(left, right) / min(left, right)
            if ratio < 3.0:
                # Fairly balanced around zero → probably already shifted
                return False

    # Default: apply the shift
    return True


# -----------------------------
# PDOS file classification
# -----------------------------

def classify_pdos_file(fname):
    """
    Classify a filename as a QE PDOS file and extract metadata.

    Valid QE PDOS patterns (projwfc.x):
    - prefix.pdos_tot[_up/down]
    - prefix.pdos_atm#N(Element)_wfc#M(orbital)[_up/down]

    Parameters
    ----------
    fname : str
        Filename in current directory.

    Returns
    -------
    info : dict or None
        If valid PDOS file:
            {
                "filename": <str>,
                "seedname": <str>,
                "kind": "tot" or "proj",
                "spin": "up", "down", or None
            }
        If not a PDOS file, returns None.
    """
    # Total DOS: prefix.pdos_tot, prefix.pdos_tot_up, prefix.pdos_tot_down
    m_tot = re.match(r"^(.+)\.pdos_tot(_(up|down))?$", fname)
    if m_tot:
        seedname = m_tot.group(1)
        spin = m_tot.group(3)  # can be 'up', 'down', or None
        return {
            "filename": fname,
            "seedname": seedname,
            "kind": "tot",
            "spin": spin,
        }

    # Projected DOS:
    # prefix.pdos_atm#1(Ti)_wfc#3(d)
    # prefix.pdos_atm#1(Ti)_wfc#3(d)_up
    m_proj = re.match(r"^(.+)\.pdos_atm#\d+\(.*?\)_wfc#\d+\(.*?\)(_up|_down)?$", fname)
    if m_proj:
        seedname = m_proj.group(1)
        spin_suffix = m_proj.group(2)  # e.g. "_up", "_down", or None
        spin = None
        if spin_suffix is not None:
            if "up" in spin_suffix:
                spin = "up"
            elif "down" in spin_suffix:
                spin = "down"
        return {
            "filename": fname,
            "seedname": seedname,
            "kind": "proj",
            "spin": spin,
        }

    # Any other file is not considered a PDOS file
    return None


# -----------------------------
# Main plotting logic
# -----------------------------

def main():
    # 1) Collect and classify all PDOS files
    entries = []
    for fname in os.listdir("."):
        info = classify_pdos_file(fname)
        if info is not None:
            entries.append(info)

    if not entries:
        print("[ERROR] No valid QE PDOS files found in current directory.")
        return

    # Sort entries for reproducible ordering
    entries = sorted(entries, key=lambda x: x["filename"])

    print("[INFO] Valid PDOS files found:")
    for e in entries:
        print("   ", e["filename"])

    # Detect whether we have spin-resolved PDOS
    has_spin = any(e["spin"] in ("up", "down") for e in entries)

    # Infer seedname from the first entry
    seedname = entries[0]["seedname"]
    print(f"[INFO] seedname detected → {seedname}")
    if has_spin:
        print("[INFO] Spin-resolved PDOS detected (up/down).")
    else:
        print("[INFO] No spin-resolved suffix found → treating as non-spin PDOS.")

    # 2) Detect Fermi energy
    ef = detect_fermi_energy_from_out()
    if ef is not None:
        print(f"[INFO] Detected Fermi energy from *.out: EF = {ef:.6f} eV")
    else:
        print("[WARN] Could not detect Fermi energy from *.out. No shift will be applied.")
        ef = None

    # 3) Decide whether to apply the Fermi shift using one representative file
    apply_fermi_shift = False
    if ef is not None:
        rep_file = entries[0]["filename"]
        try:
            E_rep, cols_rep, _ = parse_pdos_file(rep_file)
            apply_fermi_shift = should_apply_fermi_shift(E_rep, ef)
        except Exception as e:
            print(f"[WARN] Failed to parse representative PDOS file '{rep_file}': {e}")
            # If we cannot inspect, conservatively apply the shift
            apply_fermi_shift = True

        if apply_fermi_shift:
            print("[INFO] Applying E → E - EF shift.")
        else:
            print("[INFO] PDOS appears EF-centered already → no shift applied.")
    else:
        print("[INFO] EF not found → EF shift disabled.")

    # 4) Separate total and projected PDOS entries
    total_entries = [e for e in entries if e["kind"] == "tot"]
    proj_entries = [e for e in entries if e["kind"] == "proj"]

    # If we have spin-resolved PDOS, ignore non-spin entries (safety)
    if has_spin:
        total_entries = [e for e in total_entries if e["spin"] in ("up", "down")]
        proj_entries = [e for e in proj_entries if e["spin"] in ("up", "down")]

    # 5) Start plotting
    plt.figure(figsize=(8, 6))

    # x-axis
    if ef is not None and apply_fermi_shift:
        xlabel_text = "Energy (eV, shifted: E − EF)"
    else:
        xlabel_text = "Energy (eV)"
    plt.xlabel(xlabel_text)

    # y-axis depends on spin treatment
    if has_spin:
        ylabel_text = "DOS (arb. units, spin ↓ plotted negative)"
    else:
        ylabel_text = "DOS (arb. units)"
    plt.ylabel(ylabel_text)

    plt.title(f"QE PDOS overlay: {seedname}")
    plt.grid(True, alpha=0.3)

    # 5-1) Plot total DOS first (if present)
    for e in total_entries:
        fname = e["filename"]
        spin = e["spin"]  # can be 'up', 'down' or None (but here, if has_spin: up/down)

        try:
            E, cols, _ = parse_pdos_file(fname)
        except Exception as ex:
            print(f"[ERROR] Failed reading total DOS file '{fname}': {ex}")
            continue

        if ef is not None and apply_fermi_shift:
            E_plot = E - ef
        else:
            E_plot = E

        # First DOS column
        y = cols[:, 0]

        # Spin handling: down → negative
        if has_spin and spin == "down":
            y = -y
            label = extract_label_from_filename(fname) + " (↓)"
        elif has_spin and spin == "up":
            label = extract_label_from_filename(fname) + " (↑)"
        else:
            label = extract_label_from_filename(fname)

        plt.plot(E_plot, y, label=label, linewidth=2.0)

    # 5-2) Plot projected PDOS
    for e in proj_entries:
        fname = e["filename"]
        spin = e["spin"]

        try:
            E, cols, _ = parse_pdos_file(fname)
        except Exception as ex:
            print(f"[ERROR] Failed reading PDOS file '{fname}': {ex}")
            continue

        if ef is not None and apply_fermi_shift:
            E_plot = E - ef
        else:
            E_plot = E

        y = cols[:, 0]

        if has_spin and spin == "down":
            y = -y
            label = extract_label_from_filename(fname) + " (↓)"
        elif has_spin and spin == "up":
            label = extract_label_from_filename(fname) + " (↑)"
        else:
            label = extract_label_from_filename(fname)

        # alpha < 1 for readability when many curves overlay
        plt.plot(E_plot, y, label=label, linewidth=1.0, alpha=0.8)

    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()

    # 6) Output filename based on seedname
    outname = f"{seedname}_pdos_overlay.png"
    plt.savefig(outname, dpi=200)
    plt.close()

    print(f"[OK] Saved overlay plot → {outname}")


if __name__ == "__main__":
    main()

