#!/usr/bin/env bash
# Bounded serial/parallel comparison, not the full manuscript simulation.
#SBATCH --job-name=sinabro-bat-pilot
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/bat-pilot-%j.out
#SBATCH --error=logs/bat-pilot-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE_SOURCE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
PROCESSED_ROOT="${SINABRO_BAT_PROCESSED:-${PROJECT}/data/processed/bat1k_historical_v1}"
OUTPUT_ROOT="${SINABRO_BAT_PILOT_OUTPUT:-${PROJECT}/runs/bat_cosmic_pilot_${SLURM_JOB_ID}}"
N_TRAJ="${SINABRO_TEST_N_TRAJ:-100}"
MASTER_SEED="${SINABRO_SEED:-42}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); bat comparison: 1 vs $SLURM_CPUS_PER_TASK workers"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/validate_bat_pilot.py" \
  --processed-root "$PROCESSED_ROOT" \
  --out-dir "$OUTPUT_ROOT" \
  --sinabro-source "$PACKAGE_SOURCE" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --records 12 --n-traj "$N_TRAJ" --seed "$MASTER_SEED" --max-steps 10000
