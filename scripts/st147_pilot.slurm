#!/usr/bin/env bash
#SBATCH --job-name=chronoclade-st147
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=chronoclade-st147-%j.log
set -euo pipefail
# Submit from the checkout with the frozen pilot context already staged.
# Supply CHRONOCLADE_PILOT_ROOT; optionally CHRONOCLADE_PYTHON and account/partition via sbatch.
: "${CHRONOCLADE_PILOT_ROOT:?Set CHRONOCLADE_PILOT_ROOT to the staged pilot directory}"
CHRONOCLADE_PYTHON="${CHRONOCLADE_PYTHON:-python3}"
"${CHRONOCLADE_PYTHON}" scripts/st147_pilot.py --stage analyse \
  --output "${CHRONOCLADE_PILOT_ROOT}" \
  --context "${CHRONOCLADE_PILOT_ROOT}/context" \
  --threads "${SLURM_CPUS_PER_TASK:-4}" --date-randomisations 100
