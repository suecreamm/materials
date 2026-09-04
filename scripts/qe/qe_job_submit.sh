#!/bin/bash
#SBATCH --job-name=scf_pb
#SBATCH --nodelist=nodeXX          # Specific node
#SBATCH --partition=YYYY           # MPI ranks (for pw/ph/epw etc.)
#SBATCH --ntasks=64
#SBATCH --output=%x.%j.out         # Stdout file
#SBATCH --error=%x.%j.err          # Stderr file
#SBATCH --export=ALL

###############################################################################
# QE + EPW SLURM RUNNER
# (pw.x / bands.x / plotband.x / ph.x / projwfc.x / pw2wannier90.x / wannier90.x / epw.x / q2r.x / matdyn.x)
###############################################################################

# ==== 0) oneAPI environment ====
source /opt/intel/oneapi/compiler/2025.2/env/vars.sh intel64
source /opt/intel/oneapi/mkl/2025.2/env/vars.sh intel64
source /opt/intel/oneapi/mpi/2021.16/env/vars.sh intel64

# ==== 1) MPI / OpenMP settings ====
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_PIN_DOMAIN=core
export I_MPI_PIN_ORDER=compact

echo "Node list       : ${SLURM_JOB_NODELIST}"
echo "SLURM_NTASKS    : ${SLURM_NTASKS}"
NTASKS=${SLURM_NTASKS:-12}

# ==== 2) QE / EPW binaries (EDIT if your path is different) ====
QE_BIN=/home/hwang/pkgs/q-e/bin
#QE_BIN=//home/hwang/pkgs/vanila/q-e/bin

PW_BIN=${QE_BIN}/pw.x
PH_BIN=${QE_BIN}/ph.x
LAMBDA_BIN=${QE_BIN}/lambda.x
PDOS_BIN=${QE_BIN}/projwfc.x
PW2WANN_BIN=${QE_BIN}/pw2wannier90.x
EPW_BIN=${QE_BIN}/epw.x
WANN_BIN=${QE_BIN}/wannier90.x
Q2R_BIN=${QE_BIN}/q2r.x
MATDYN_BIN=${QE_BIN}/matdyn.x

BANDS_BIN=${QE_BIN}/bands.x
PLOTBAND_BIN=${QE_BIN}/plotband.x

# ==== 3) Input lists (FILL THESE ARRAYS) =======================

# (A) Ground-state / standard pw.x runs (SCF/NSCF coarse meshes etc.)
PW_INPUTS=(
1scf.in
#3nscf.in
)

BAND_PW_INPUTS=(
#99band.1pw.in
)

BANDS_INPUTS=(
#99band.2pp.in
)

# (D) Optional: plotband.x inputs (makes plottable data)
PLOTBAND_INPUTS=(
)

# phonons
PH_INPUTS=(
#2ph.in
)

PDOS_INPUTS=(
#1scf.2pdos.in
)

PW2WAN_INPUTS=(
)

WANN_SEEDS=(
)

EPW_INPUTS=(
#4.1epw1.in
#4.2epw2.in
#4.3epw3.in
)

Q2R_INPUTS=(
#2ph.1q2r.in
)

MATDYN_INPUTS=(
#2ph.2matdyn_dos.in
#2ph.2matdyn_band.in
)

echo "======================================"
echo "PW_INPUTS        : ${PW_INPUTS[*]}"
echo "BAND_PW_INPUTS   : ${BAND_PW_INPUTS[*]}"
echo "BANDS_INPUTS     : ${BANDS_INPUTS[*]}"
echo "PLOTBAND_INPUTS  : ${PLOTBAND_INPUTS[*]}"
echo "PH_INPUTS        : ${PH_INPUTS[*]}"
echo "PDOS_INPUTS      : ${PDOS_INPUTS[*]}"
echo "PW2WAN_INPUTS    : ${PW2WAN_INPUTS[*]}"
echo "WANN_SEEDS       : ${WANN_SEEDS[*]}"
echo "EPW_INPUTS        : ${EPW_INPUTS[*]}"
echo "Q2R_INPUTS        : ${Q2R_INPUTS[*]}"
echo "MATDYN_INPUTS     : ${MATDYN_INPUTS[*]}"
echo "======================================"

# ==== 4) Run pw.x (SCF/NSCF etc.) ===============================
if ((${#PW_INPUTS[@]})); then
  echo ">>> Running pw.x jobs (standard)"
  for inpf in "${PW_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PW] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PW] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${PW_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No PW_INPUTS; skipping pw.x (standard) block"
fi

###############################################################################
# ==== 5) Band structure workflow (BEFORE ph.x) ===============================
###############################################################################

# 5-1) pw.x on band-path (calculation='bands' OR 'nscf' with band-path K_POINTS)
if ((${#BAND_PW_INPUTS[@]})); then
  echo ">>> Running pw.x band-path jobs"
  for inpf in "${BAND_PW_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PW-BAND] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PW-BAND] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${PW_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No BAND_PW_INPUTS; skipping pw.x band-path block"
fi

# 5-2) bands.x (single-core) - extracts band energies for plotting/analysis
if ((${#BANDS_INPUTS[@]})); then
  echo ">>> Running bands.x jobs (single-core)"
  for inpf in "${BANDS_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[BANDS] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[BANDS] ${inpf} -> ${stem}.out"
    srun -n 1 "${BANDS_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No BANDS_INPUTS; skipping bands.x block"
fi

# 5-3) plotband.x (optional, single-core)
if ((${#PLOTBAND_INPUTS[@]})); then
  echo ">>> Running plotband.x jobs (single-core)"
  for inpf in "${PLOTBAND_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PLOTBAND] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PLOTBAND] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n ""${NTASKS} "${PLOTBAND_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No PLOTBAND_INPUTS; skipping plotband.x block"
fi

# ==== 6) Run ph.x ==============================================
if ((${#PH_INPUTS[@]})); then
  echo ">>> Running ph.x jobs"
  for inpf in "${PH_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PH] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PH] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${PH_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No PH_INPUTS; skipping ph.x block"
fi

# ==== 7) Run projwfc.x (PDOS) =================================
if ((${#PDOS_INPUTS[@]})); then
  echo ">>> Running projwfc.x (PDOS) jobs"
  for inpf in "${PDOS_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PDOS] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PDOS] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n 4 "${PDOS_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No PDOS_INPUTS; skipping projwfc.x block"
fi

# ==== 8) Run pw2wannier90.x (MPI OK) ============================
if ((${#PW2WAN_INPUTS[@]})); then
  echo ">>> Running pw2wannier90.x jobs"
  for inpf in "${PW2WAN_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[PW2WANN] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[PW2WANN] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${PW2WANN_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No PW2WAN_INPUTS; skipping pw2wannier90.x block"
fi

# ==== 9) Run wannier90.x (-pp then run) =========================
if ((${#WANN_SEEDS[@]})); then
  echo ">>> Running wannier90.x jobs (single-core: -pp then main)"
  for seed in "${WANN_SEEDS[@]}"; do
    if [[ ! -f "${seed}.win" ]]; then
      echo "[WANN] WARNING: '${seed}.win' not found, skipping"
      continue
    fi
    echo "[WANN-PP] ${seed} (-pp)"
    srun -n 1 "${WANN_BIN}" -pp "${seed}" > "${seed}.pp.out"

    echo "[WANN] ${seed} (main run)"
    srun -n 1 "${WANN_BIN}" "${seed}" > "${seed}.wout"
  done
else
  echo ">>> No WANN_SEEDS; skipping wannier90.x block"
fi

# ==== 9.5) Link ph outputs before EPW ===========================
if ((${#EPW_INPUTS[@]})); then
  echo ">>> Linking ph outputs for EPW"

  LINK_SCRIPT="$(dirname "$0")/../link_ph_outputs.py"   # 4epw/../link_ph_outputs.py
  PREFIX="TiS2"                                          # ★ 여기만 수정
  PH_ROOT=".."
  DVSCF_DIR="./tmp/_ph0"

  if [[ ! -f "${LINK_SCRIPT}" ]]; then
    echo "[LINK] ERROR: link script not found at ${LINK_SCRIPT}"
    exit 1
  fi

  python3 "${LINK_SCRIPT}" "${PREFIX}" \
    --ph_root "${PH_ROOT}" \
    --dvscf_dir "${DVSCF_DIR}"

  if [[ $? -ne 0 ]]; then
    echo "[LINK] ERROR: linking failed, aborting EPW"
    exit 1
  fi
fi

# ==== 10) Run epw.x =============================================
if ((${#EPW_INPUTS[@]})); then
  echo ">>> Running epw.x jobs"
  for inpf in "${EPW_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[EPW] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[EPW] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${EPW_BIN}" -npool "${NTASKS}" -nimage 1 -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No EPW_INPUTS; skipping epw.x block"
fi

# ==== 11) Run q2r.x (single-core) ===============================
if ((${#Q2R_INPUTS[@]})); then
  echo ">>> Running q2r.x jobs (single-core)"
  for inpf in "${Q2R_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[Q2R] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[Q2R] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${Q2R_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No Q2R_INPUTS; skipping q2r.x block"
fi

# ==== 12) Run matdyn.x (single-core) ===========================
if ((${#MATDYN_INPUTS[@]})); then
  echo ">>> Running matdyn.x jobs (single-core)"
  for inpf in "${MATDYN_INPUTS[@]}"; do
    if [[ ! -f "$inpf" ]]; then
      echo "[MATDYN] WARNING: input file '$inpf' not found, skipping"
      continue
    fi
    stem="${inpf%.in}"
    echo "[MATDYN] ${inpf} -> ${stem}.out"
    srun --mpi=pmi2 -n "${NTASKS}" "${MATDYN_BIN}" -in "${inpf}" > "${stem}.out"
  done
else
  echo ">>> No MATDYN_INPUTS; skipping matdyn.x block"
fi

echo "=== All requested jobs finished ==="


