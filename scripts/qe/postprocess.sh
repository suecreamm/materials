"""
Purpose:
Slurm batch script for running EPW post-processing and plotting
using a Python workflow.

What this script does:
  1) Runs on the 'master' partition with a single task.
  2) Activates a conda environment ('pymatgen') in a batch-safe way.
  3) Executes the Python post-processing script:
      postprocess_epw.py <PREFIX>
Usage:
  sbatch postprocess.sh

Requirements:
  - Conda must be available on the compute node.
  - A conda environment named 'pymatgen' must exist.
  - The file 'postprocess_epw.py' must be in the same directory
    (or in the PYTHONPATH).
Notes:
  - This script is intended to be used together with postprocess_epw.py.
  - PREFIX (e.g., graphene) should match the EPW/QE calculation prefix.
"""

#!/bin/bash
#SBATCH --job-name=EPWpostprocess
#SBATCH --nodelist=master          # select the node
#SBATCH --partition=master         # select the node
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --export=ALL

set -euo pipefail

PREFIX="graphene"

echo "=== Postprocess start ==="
echo "Host: $(hostname)"
echo "PWD : $(pwd)"
date

# --- Batch-safe conda activate ---
if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda not found in PATH on this node."
  exit 2
fi

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate pymatgen                          # conda environment

echo "[info] python: $(which python)"
python -c "import sys; print('[info] sys.executable=', sys.executable)"
python -c "import pymatgen; print('[info] pymatgen OK')"

# --- Run your plotting / postprocess ---
python postprocess_epw.py "${PREFIX}"            # code to run

echo "=== Postprocess done ==="
date

