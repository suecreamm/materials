#!/bin/bash
#SBATCH --job-name=epw1
#SBATCH --partition=7501
#SBATCH --ntasks=64
#SBATCH --output=qe_full.out
#SBATCH --error=qe_full.err
#SBATCH --export=ALL
#SBATCH --hint=nomultithread

###############################################################################
# EPW runner from current directory (e.g. retry/4.1epw, retry/4.2epw, ...)
#
# Required layout:
#   current dir : .../pb/4.1epw
#   parent dir  : .../pb
#   source tmp  : ../tmp
#   source dyn  : ../pb.dyn*
#
# What this script does:
#   1) Initializes oneAPI environment in the same style as the known-working job
#   2) Creates ./tmp if missing
#   3) Relative-links ../tmp/* into ./tmp/
#   4) Relative-links ../pb.dyn* into current directory
#   5) Creates/refreshes EPW-required links in ./tmp/_ph0
#   6) Runs epw.x for input(s) listed in EPW_INPUTS
###############################################################################

ulimit -s unlimited
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Keep the environment initialization consistent with the known-working script.
if [[ -z "${ONEAPI_ROOT}" && -f /opt/intel/oneapi/setvars.sh ]]; then
  source /opt/intel/oneapi/setvars.sh
fi

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=core
export I_MPI_PIN_ORDER=compact

echo "Running with NTASKS=${SLURM_NTASKS:-64}"

QE_BIN_DIR=/home/hwang/pkgs/q-e/bin
EPW_BIN="${QE_BIN_DIR}/epw.x"
NTASKS=${SLURM_NTASKS:-64}

PREFIX="pb"

# ------------------------------------------------------------------
# Choose EPW inputs to run in THIS directory
# ------------------------------------------------------------------
EPW_INPUTS=(
#4.1epw.in
#4.2epw.in
#4.3epw.in
#4.4epw.in
4.5epw.in
)

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
WORKDIR="$(pwd)"
ROOTDIR="$(cd .. && pwd)"
SRC_TMP="${ROOTDIR}/tmp"
DST_TMP="${WORKDIR}/tmp"

echo "======================================"
echo "WORKDIR    : ${WORKDIR}"
echo "ROOTDIR    : ${ROOTDIR}"
echo "SRC_TMP    : ${SRC_TMP}"
echo "DST_TMP    : ${DST_TMP}"
echo "PREFIX     : ${PREFIX}"
echo "EPW_BIN    : ${EPW_BIN}"
echo "NTASKS     : ${NTASKS}"
echo "NPOOL      : ${NPOOL}"
echo "EPW_INPUTS : ${EPW_INPUTS[*]}"
echo "======================================"

# ------------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------------
if [[ ! -x "${EPW_BIN}" ]]; then
  echo "[ERROR] EPW binary not found or not executable: ${EPW_BIN}"
  exit 1
fi

if [[ ! -d "${SRC_TMP}" ]]; then
  echo "[ERROR] Source tmp directory does not exist: ${SRC_TMP}"
  exit 1
fi

# ------------------------------------------------------------------
# Utility: relative symlink
# ------------------------------------------------------------------
relative_link() {
  local src="$1"
  local dst="$2"

  if [[ ! -e "${src}" && ! -L "${src}" ]]; then
    echo "[ERROR] Source does not exist: ${src}"
    return 1
  fi

  mkdir -p "$(dirname "${dst}")"

  if [[ -e "${dst}" && ! -L "${dst}" ]]; then
    echo "[SKIP] Regular existing path kept: ${dst}"
    return 0
  fi

  if [[ -L "${dst}" ]]; then
    rm -f "${dst}"
  fi

  local rel_src
  rel_src="$(realpath --relative-to="$(dirname "${dst}")" "${src}")"
  ln -s "${rel_src}" "${dst}"
  echo "[LINK] ${dst} -> ${rel_src}"
}

# ------------------------------------------------------------------
# 1) Ensure ./tmp exists
# ------------------------------------------------------------------
mkdir -p "${DST_TMP}"

# ------------------------------------------------------------------
# 2) Relative-link ../tmp/* into ./tmp/
# ------------------------------------------------------------------
echo ">>> Linking ${SRC_TMP} contents into ${DST_TMP}"

shopt -s nullglob dotglob
for src in "${SRC_TMP}"/*; do
  base="$(basename "${src}")"
  dst="${DST_TMP}/${base}"
  relative_link "${src}" "${dst}" || exit 1
done
shopt -u nullglob dotglob

# Ensure _ph0 exists as a real directory if absent
mkdir -p "${DST_TMP}/_ph0"

# ------------------------------------------------------------------
# 3) Relative-link ../PREFIX.dyn* into current directory
# ------------------------------------------------------------------
echo ">>> Linking ${ROOTDIR}/${PREFIX}.dyn* into current directory"

shopt -s nullglob
dyn_candidates=( "${ROOTDIR}/${PREFIX}".dyn* )
shopt -u nullglob

if [[ ${#dyn_candidates[@]} -eq 0 ]]; then
  echo "[ERROR] No dyn files found in ${ROOTDIR}"
  exit 1
fi

for src in "${dyn_candidates[@]}"; do
  base="$(basename "${src}")"
  dst="${WORKDIR}/${base}"
  relative_link "${src}" "${dst}" || exit 1
done

# ------------------------------------------------------------------
# 4) Build EPW-required links inside ./tmp/_ph0
# ------------------------------------------------------------------
echo ">>> Preparing dyn_q / dvscf_q links under ./tmp/_ph0"

python3 - "${PREFIX}" "${DST_TMP}/_ph0" <<'PYCODE'
from __future__ import annotations

from pathlib import Path
import os
import re
import sys


def safe_symlink(src: Path, dst: Path) -> None:
    """
    Replace dst with a symlink to src.
    Keep regular files untouched.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() and not dst.is_symlink():
        print(f"[skip] regular file kept: {dst}")
        return

    if dst.is_symlink() or dst.exists():
        dst.unlink()

    rel = os.path.relpath(src.resolve(), start=dst.parent.resolve())
    os.symlink(rel, dst)
    print(f"[link] {dst} -> {rel}")


def first_existing(candidates: list[Path]) -> Path | None:
    """
    Return the first existing regular file from candidates.
    """
    for f in candidates:
        if f.is_file():
            return f
    return None


def main(prefix: str, ph0_dir: Path) -> int:
    cwd = Path.cwd()
    ph0_dir = ph0_dir.resolve()

    print(f"[info] cwd    : {cwd}")
    print(f"[info] ph0_dir: {ph0_dir}")

    # --------------------------------------------------------------
    # Link PREFIX.dynN -> ./tmp/_ph0/PREFIX.dyn_qN
    # --------------------------------------------------------------
    dyn_pat = re.compile(rf"^{re.escape(prefix)}\.dyn(\d+)$")
    dyn_hits = 0

    for f in sorted(cwd.glob(f"{prefix}.dyn*")):
        m = dyn_pat.match(f.name)
        if not m:
            continue

        iq = m.group(1)
        safe_symlink(f, ph0_dir / f"{prefix}.dyn_q{iq}")
        dyn_hits += 1

    if dyn_hits == 0:
        print("[error] no usable dyn files found in current directory")
        return 1

    # --------------------------------------------------------------
    # Search dvscf source files
    #
    # Standard QE ph.x layout with:
    #
    #   prefix   = 'pb'
    #   fildvscf = 'dvscf'
    #
    # is normally:
    #
    #   _ph0/pb.dvscf1
    #   _ph0/pb.q_2/pb.dvscf1
    #   _ph0/pb.q_3/pb.dvscf1
    #   ...
    #
    # Some previous calculations may instead contain:
    #
    #   pb.pb.dvscf1
    #
    # so that form is accepted as a fallback.
    #
    # IMPORTANT:
    # For q >= 2, the q-index comes from the directory name
    # "pb.q_N", NOT from the trailing "1" in "pb.dvscf1".
    # --------------------------------------------------------------

    hits: dict[int, Path] = {}

    # --------------------------------------------------------------
    # q1 / Gamma
    #
    # Prefer the standard file produced by fildvscf='dvscf'.
    # Keep the old pb.pb.dvscf1 naming as a fallback.
    # --------------------------------------------------------------
    q1_candidates = [
        ph0_dir / f"{prefix}.dvscf1",
        ph0_dir / f"{prefix}.{prefix}.dvscf1",
    ]

    q1_src = first_existing(q1_candidates)

    if q1_src is not None:
        hits[1] = q1_src
        print(f"[info] q1 dvscf source: {q1_src}")
    else:
        print("[warn] q1 dvscf source not found directly under _ph0")

    # --------------------------------------------------------------
    # q2, q3, ...
    #
    # QE directories:
    #   pb.q_2/
    #   pb.q_3/
    #   ...
    # --------------------------------------------------------------
    qdir_pat = re.compile(
        rf"^{re.escape(prefix)}\.q_(\d+)$"
    )

    for qdir in sorted(ph0_dir.glob(f"{prefix}.q_*")):
        if not qdir.is_dir():
            continue

        m = qdir_pat.match(qdir.name)
        if not m:
            continue

        iq = int(m.group(1))

        # q1 is normally handled from _ph0 directly.
        # If pb.q_1 exists, only use it as a fallback when no
        # direct q1 dvscf file was found.
        if iq == 1 and 1 in hits:
            continue

        candidates = [
            # Standard form for:
            # prefix='pb', fildvscf='dvscf'
            qdir / f"{prefix}.dvscf1",

            # Older / exceptional form for:
            # fildvscf='pb.dvscf'
            qdir / f"{prefix}.{prefix}.dvscf1",
        ]

        src = first_existing(candidates)

        if src is None:
            print(
                f"[warn] no dvscf source found for q{iq} in {qdir}"
            )
            continue

        hits[iq] = src
        print(f"[info] q{iq} dvscf source: {src}")

    if not hits:
        print("[error] no dvscf source files found")
        return 2

    print(
        f"[info] found dvscf q-points: {sorted(hits.keys())}"
    )

    # --------------------------------------------------------------
    # Build EPW-required links
    #
    # Keep BOTH historical link forms because the existing runner
    # already created them and other workflows may rely on them.
    # --------------------------------------------------------------
    for iq in sorted(hits):
        src = hits[iq]

        safe_symlink(
            src,
            ph0_dir / f"{prefix}.dvscf{iq}_1"
        )

        safe_symlink(
            src,
            ph0_dir / f"{prefix}.dvscf_q{iq}"
        )

    print("[done] EPW link preparation finished")
    return 0


if __name__ == "__main__":
    prefix = sys.argv[1]
    ph0_dir = Path(sys.argv[2])
    sys.exit(main(prefix, ph0_dir))
PYCODE

rc=$?
if [[ $rc -ne 0 ]]; then
  echo "[ERROR] EPW pre-link step failed with exit code ${rc}"
  exit $rc
fi

# ------------------------------------------------------------------
# 5) Debug listing
# ------------------------------------------------------------------
echo ">>> Debug listing: current directory"
ls -lh

echo ">>> Debug listing: ./tmp"
ls -lh ./tmp | sed -n '1,120p'

echo ">>> Debug listing: ./tmp/_ph0"
ls -lh ./tmp/_ph0 | sed -n '1,220p' || true

# ------------------------------------------------------------------
# 6) Run EPW
# ------------------------------------------------------------------
if ((${#EPW_INPUTS[@]})); then
  echo ">>> Running epw.x jobs"
  for inpf in "${EPW_INPUTS[@]}"; do
    if [[ ! -f "${inpf}" ]]; then
      echo "[EPW] WARNING: input file '${inpf}' not found, skipping"
      continue
    fi

    stem="${inpf%.in}"
    echo "[EPW] ${inpf} -> ${stem}.out"

    # Use -npool NPOOL, not -npool NTASKS unless you intentionally want that.
    srun --mpi=pmi2 -n "${NTASKS}" "${EPW_BIN}" -npool "${NTASKS}" -nimage 1 -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No EPW_INPUTS; skipping epw.x block"
fi

echo "=== All requested jobs finished ==="
