#!/usr/bin/env bash
# Read every saved path from all three datasets, without resimulation.
#SBATCH --job-name=sinabro-apobec-audit
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/apobec-audit-%j.out
#SBATCH --error=logs/apobec-audit-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
INPUTS="${SINABRO_APOBEC_INPUTS:-${PROJECT}/data/processed/experimental_apobec_inputs_v1}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
SEED="${SINABRO_SEED:-42}"
RUN_ROOT="${SINABRO_APOBEC_RUN:-${PROJECT}/runs/hpc_experimental_apobec_n${N_TRAJ}_seed${SEED}}"
OUTPUT="${SINABRO_APOBEC_AUDIT_OUTPUT:-${PROJECT}/runs/experimental_apobec_audit_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); auditing experimental APOBEC results: $RUN_ROOT"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/validate_campaign.py" \
  --inputs-dir "$INPUTS" --run-root "$RUN_ROOT" --out-root "$OUTPUT" \
  --n-traj "$N_TRAJ" --workers "$SLURM_CPUS_PER_TASK"
