#!/usr/bin/env bash
# Descriptive ranks and distribution figures from an audited table export.
#SBATCH --job-name=sinabro-distributions
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/distributions-%j.out
#SBATCH --error=logs/distributions-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
: "${SINABRO_DISTRIBUTIONS_INPUT:?Export the completed aggregation_JOBID directory before sbatch}"

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
OUTPUT_DIR="${SINABRO_DISTRIBUTIONS_OUTPUT:-${PROJECT}/results/signature_distributions_${SLURM_JOB_ID}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); analyzing exported tables: $SINABRO_DISTRIBUTIONS_INPUT"
"$PYTHON_BIN" "$ANALYSIS_DIR/03_mean_std_analysis/analyze_signature_distributions.py" \
  --aggregation-dir "$SINABRO_DISTRIBUTIONS_INPUT" \
  --out-dir "$OUTPUT_DIR" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --expected-records 73214 \
  --expected-signatures 79
