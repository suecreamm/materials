#!/usr/bin/env python3
"""
Plot EPW electron self-energy along a band path.

Input file format (elecselfen = .true., 4 columns):
    ik   ibnd   E(ibnd) [eV, relative to E_F]   Im(Sigma) [meV]

Note: Re(Sigma) is NOT in this file. It is printed in the EPW stdout
(4.5epw.out). See parse_resigma_from_out() at the bottom.

Usage:
    python3 plot_elself.py
    python3 plot_elself.py linewidth.elself.0.075K pb_band.kpt
"""

import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HBAR_MEV_FS = 658.2119   # hbar in meV*fs

# ------------------------------------------------------------------
# 1. read the self-energy file
# ------------------------------------------------------------------
if len(sys.argv) > 1:
    self_file = sys.argv[1]
else:
    cand = sorted(glob.glob("linewidth.elself*"))
    if not cand:
        sys.exit("linewidth.elself.* not found")
    self_file = cand[0]

kpt_file = sys.argv[2] if len(sys.argv) > 2 else "pb_band.kpt"

data = np.loadtxt(self_file, comments="#")
if data.shape[1] != 4:
    sys.exit("expected 4 columns (ik ibnd E ImSigma), got %d" % data.shape[1])

ik    = data[:, 0].astype(int)
ibnd  = data[:, 1].astype(int)
ener  = data[:, 2]              # eV, relative to E_F
imsig = data[:, 3]              # meV

nk    = ik.max()
bands = np.unique(ibnd)

# Im(Sigma) == 0 means "not computed" (k point outside fsthick),
# not "zero linewidth". Mask it so it is never drawn as a real point.
computed = imsig > 0.0
gamma = np.where(computed, 2.0 * imsig, np.nan)                    # meV
tau = np.full_like(gamma, np.nan)                                  # fs
tau[computed] = HBAR_MEV_FS / gamma[computed]

print("file        : %s" % self_file)
print("k points    : %d,  bands : %s" % (nk, list(bands)))
print("computed    : %d / %d states inside fsthick"
      % (computed.sum(), len(imsig)))

# ------------------------------------------------------------------
# 2. high symmetry ticks from the filkf file
# ------------------------------------------------------------------
HIGHSYM = {
    (0.000, 0.000, 0.000): "G",
    (0.500, 0.000, 0.500): "X",
    (0.500, 0.250, 0.750): "W",
    (0.500, 0.500, 0.500): "L",
    (0.375, 0.375, 0.750): "K",
    (0.625, 0.250, 0.625): "U",
}

ticks, labels = [], []
try:
    with open(kpt_file) as f:
        lines = [l for l in f if l.strip()]
    kpts = np.array([[float(x) for x in l.split()[:3]]
                     for l in lines[1:] if len(l.split()) >= 3])
    for i, k in enumerate(kpts):
        key = tuple(np.round(np.abs(k), 3))
        if key in HIGHSYM and (not ticks or ticks[-1] != i + 1):
            ticks.append(i + 1)
            labels.append(HIGHSYM[key])
    print("ticks       : %s" % list(zip(labels, ticks)))
except (IOError, ValueError):
    print("ticks       : k path file not usable, using indices")


def decorate(ax):
    for t in ticks:
        ax.axvline(t, color="0.75", lw=0.7, zorder=0)
    if ticks:
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    ax.set_xlim(1, nk)


# ------------------------------------------------------------------
# 3. band structure, marker size and colour = linewidth
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5.5))

for b in bands:
    m = ibnd == b
    o = np.argsort(ik[m])
    ax.plot(ik[m][o], ener[m][o], "-", color="0.8", lw=1.0, zorder=1)

m = computed
sc = ax.scatter(ik[m], ener[m], c=gamma[m], s=8 + 3.0 * gamma[m],
                cmap="plasma", edgecolors="none", zorder=3)
cb = fig.colorbar(sc, ax=ax)
cb.set_label("Linewidth  2 Im(Sigma)  (meV)")

ax.axhline(0.0, color="k", lw=0.9, ls="--", zorder=2)
ax.axhline(1.0, color="tab:blue", lw=0.7, ls=":", zorder=2)
ax.axhline(-1.0, color="tab:blue", lw=0.7, ls=":", zorder=2)
ax.text(2.0, 1.05, "fsthick window", color="tab:blue", fontsize=8)

ax.set_ylabel("E - E_F  (eV)")
ax.set_title("Pb: bands with electron-phonon linewidth  (T = 0.075 K)")
decorate(ax)
fig.tight_layout()
fig.savefig("elself_bands.png", dpi=200)
print("wrote elself_bands.png")

# ------------------------------------------------------------------
# 4. linewidth and lifetime along the path (masked, no fake zeros)
# ------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

for b in bands:
    m = ibnd == b
    o = np.argsort(ik[m])
    axes[0].plot(ik[m][o], gamma[m][o], ".-", ms=4, lw=1.0,
                 label="band %d" % b)
    axes[1].plot(ik[m][o], tau[m][o], ".-", ms=4, lw=1.0,
                 label="band %d" % b)

axes[0].set_ylabel("2 Im(Sigma)  (meV)")
axes[0].legend(fontsize=8)
axes[1].set_ylabel("lifetime  hbar / 2 Im(Sigma)  (fs)")
axes[1].set_yscale("log")
axes[1].set_xlabel("k path")
for a in axes:
    decorate(a)
fig.suptitle("Gaps: k points outside fsthick, not computed", fontsize=9)
fig.tight_layout()
fig.savefig("elself_path.png", dpi=200)
print("wrote elself_path.png")

# ------------------------------------------------------------------
# 5. linewidth vs energy
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
for b in bands:
    m = computed & (ibnd == b)
    ax.plot(ener[m], gamma[m], "o", ms=4, alpha=0.75, label="band %d" % b)
ax.axvline(0.0, color="k", lw=0.9, ls="--")
ax.set_xlabel("E - E_F  (eV)")
ax.set_ylabel("2 Im(Sigma)  (meV)")
ax.set_title("Linewidth vs energy")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("elself_vs_energy.png", dpi=200)
print("wrote elself_vs_energy.png")

# ------------------------------------------------------------------
# 6. summary near E_F
# ------------------------------------------------------------------
near = computed & (np.abs(ener) < 0.05)
if near.any():
    print("\nstates within 50 meV of E_F : %d" % near.sum())
    print("  mean linewidth : %8.3f meV" % gamma[near].mean())
    print("  max  linewidth : %8.3f meV" % gamma[near].max())
    print("  mean lifetime  : %8.2f fs" % tau[near].mean())


# ------------------------------------------------------------------
# optional: Re(Sigma) from the EPW stdout
# ------------------------------------------------------------------
def parse_resigma_from_out(out_file):
    """EPW prints Re[Sigma], Im[Sigma], Z and lambda to stdout.
    The line layout is version dependent, so pull numbers by keyword."""
    import re
    rows = []
    cur_ik = cur_ib = None
    with open(out_file) as f:
        for line in f:
            m = re.search(r"ik\s*=\s*(\d+).*ibnd\s*=\s*(\d+)", line)
            if m:
                cur_ik, cur_ib = int(m.group(1)), int(m.group(2))
                continue
            if "Re[Sigma]" in line and cur_ik is not None:
                nums = re.findall(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", line)
                nums = [float(x.replace("D", "E").replace("d", "e"))
                        for x in nums]
                if len(nums) >= 3:
                    rows.append([cur_ik, cur_ib] + nums[:4])
    return np.array(rows) if rows else None
