#!/usr/bin/env bash
# Export already audited paths; no trajectory simulation or plotting.
#SBATCH --job-name=sinabro-apobec-export
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/apobec-export-%j.out
#SBATCH --error=logs/apobec-export-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
: "${SINABRO_APOBEC_AUDIT:?Set the passing experimental_apobec_audit_JOBID directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
INPUTS="${SINABRO_APOBEC_INPUTS:-${PROJECT}/data/processed/experimental_apobec_inputs_v1}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
SEED="${SINABRO_SEED:-42}"
RUN_ROOT="${SINABRO_APOBEC_RUN:-${PROJECT}/runs/hpc_experimental_apobec_n${N_TRAJ}_seed${SEED}}"
OUTPUT="${SINABRO_APOBEC_EXPORT_OUTPUT:-${PROJECT}/runs/experimental_apobec_export_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); exporting audited experimental APOBEC results: $RUN_ROOT"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/export_results.py" \
  --inputs-dir "$INPUTS" --run-root "$RUN_ROOT" --audit-root "$SINABRO_APOBEC_AUDIT" \
  --out-root "$OUTPUT" --n-traj "$N_TRAJ" --workers "$SLURM_CPUS_PER_TASK"
