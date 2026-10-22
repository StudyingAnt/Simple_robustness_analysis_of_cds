#!/usr/bin/env bash
# Original Figures S1/S2: 79 profiles x 10,000 random sequences, one event each.
#SBATCH --job-name=sinabro-paper-sampler
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/paper-sampler-%j.out
#SBATCH --error=logs/paper-sampler-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
OUTPUT="${SINABRO_SAMPLER_OUTPUT:-${PROJECT}/results/signature_sampler_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); original S1/S2 sampler validation"
"$PYTHON_BIN" "$ANALYSIS_DIR/07_manuscript_validation/validate_signature_sampler.py" \
  --sinabro-source "$PACKAGE" --out-dir "$OUTPUT" --workers "$SLURM_CPUS_PER_TASK" \
  --n-sequences 10000 --length 2200 --gc 0.4 --seed 42
