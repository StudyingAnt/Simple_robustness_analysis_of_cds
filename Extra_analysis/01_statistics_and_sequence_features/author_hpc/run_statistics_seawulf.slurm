#!/usr/bin/env bash
# Author SeaWulf configuration; Python remains scheduler-independent.
#SBATCH --job-name=sinabro-extra-stats
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/extra-statistics-%j.out
#SBATCH --error=logs/extra-statistics-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch from the project directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
HUMAN="${SINABRO_EXTRA_HUMAN:-${PROJECT}/runs/aggregation_2305071}"
BATS="${SINABRO_EXTRA_BATS:-${PROJECT}/runs/bat_cosmic_export_2328086}"
APOBEC="${SINABRO_EXTRA_APOBEC:-${PROJECT}/runs/experimental_apobec_export_2338741}"
OUTPUT="${SINABRO_EXTRA_OUTPUT:-${PROJECT}/results/extra_analysis/statistics_${SLURM_JOB_ID}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); additional paired statistics and sequence features"
"$PYTHON_BIN" "$ANALYSIS_DIR/Extra_analysis/01_statistics_and_sequence_features/run_analysis.py" \
  --human-export "$HUMAN" --bat-export "$BATS" --apobec-export "$APOBEC" \
  --out-dir "$OUTPUT" --n-traj 1000
