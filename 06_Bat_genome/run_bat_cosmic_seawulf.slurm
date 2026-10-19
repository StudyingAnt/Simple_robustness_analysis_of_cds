#!/usr/bin/env bash
# 100 shards per species, interleaved; at most 8 single-node jobs at once.
#SBATCH --job-name=sinabro-bat-cosmic
#SBATCH --partition=short-96core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --array=0-199%8
#SBATCH --time=04:00:00
#SBATCH --output=logs/bat-cosmic-%A_%a.out
#SBATCH --error=logs/bat-cosmic-%A_%a.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_ARRAY_TASK_ID:?Requires a Slurm array}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
: "${SINABRO_BAT_PILOT:?Export the passing bat_cosmic_pilot_JOBID directory before sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE_SOURCE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
PROCESSED_ROOT="${SINABRO_BAT_PROCESSED:-${PROJECT}/data/processed/bat1k_historical_v1}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
MASTER_SEED="${SINABRO_SEED:-42}"
MAX_STEPS="${SINABRO_MAX_STEPS:-10000}"
OUTPUT_ROOT="${SINABRO_BAT_RUN:-${PROJECT}/runs/hpc_bat_cosmic_n${N_TRAJ}_seed${MASTER_SEED}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID on $(hostname); bat COSMIC production: $OUTPUT_ROOT"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/run_bat_cosmic.py" \
  --processed-root "$PROCESSED_ROOT" --pilot-dir "$SINABRO_BAT_PILOT" \
  --out-root "$OUTPUT_ROOT" --array-index "$SLURM_ARRAY_TASK_ID" --shards-per-species 100 \
  --sinabro-source "$PACKAGE_SOURCE" --workers "$SLURM_CPUS_PER_TASK" \
  --n-traj "$N_TRAJ" --seed "$MASTER_SEED" --max-steps "$MAX_STEPS"
