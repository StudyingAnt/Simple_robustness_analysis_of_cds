#!/usr/bin/env bash
# Submit from the HPC project root; create logs/ before submission.
#SBATCH --job-name=sinabro-rums
#SBATCH --partition=short-96core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --array=0-99%8
#SBATCH --time=04:00:00
#SBATCH --output=logs/rums-%A_%a.out
#SBATCH --error=logs/rums-%A_%a.err

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
: "${SLURM_ARRAY_TASK_ID:?This script requires a Slurm array}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"

# The submission directory is the project root unless explicitly overridden.
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE_SOURCE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
INPUT_FASTA="${SINABRO_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
MASTER_SEED="${SINABRO_SEED:-42}"
MAX_STEPS="${SINABRO_MAX_STEPS:-10000}"
OUTPUT_ROOT="${SINABRO_OUTPUT:-${PROJECT}/runs/hpc_gencode_v40_n${N_TRAJ}_seed${MASTER_SEED}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID on $(hostname)"
echo "Input: $INPUT_FASTA; output: $OUTPUT_ROOT; paths per combination: $N_TRAJ"

"$PYTHON_BIN" "$ANALYSIS_DIR/02_run_sinabro/run_rums.py" \
  --input "$INPUT_FASTA" \
  --out-dir "$OUTPUT_ROOT" \
  --sinabro-source "$PACKAGE_SOURCE" \
  --n-traj "$N_TRAJ" \
  --seed "$MASTER_SEED" \
  --max-steps "$MAX_STEPS" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --n-shards 100 \
  --expected-records 73214 \
  --shard-index "$SLURM_ARRAY_TASK_ID"
