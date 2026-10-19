#!/usr/bin/env bash
# Export audited transcript-level tables, without any new simulation.
#SBATCH --job-name=sinabro-bat-export
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/bat-export-%j.out
#SBATCH --error=logs/bat-export-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
: "${SINABRO_BAT_AUDIT:?Export the passing bat_cosmic_audit_JOBID directory before sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PROCESSED_ROOT="${SINABRO_BAT_PROCESSED:-${PROJECT}/data/processed/bat1k_historical_v1}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
MASTER_SEED="${SINABRO_SEED:-42}"
RUN_ROOT="${SINABRO_BAT_RUN:-${PROJECT}/runs/hpc_bat_cosmic_n${N_TRAJ}_seed${MASTER_SEED}}"
OUTPUT_ROOT="${SINABRO_BAT_EXPORT_OUTPUT:-${PROJECT}/runs/bat_cosmic_export_${SLURM_JOB_ID}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); exporting both bats: $RUN_ROOT"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/process_bat_results.py" \
  --stage export --processed-root "$PROCESSED_ROOT" --run-root "$RUN_ROOT" \
  --audit-root "$SINABRO_BAT_AUDIT" --out-root "$OUTPUT_ROOT" \
  --n-traj "$N_TRAJ" --workers "$SLURM_CPUS_PER_TASK"
