#!/usr/bin/env bash
# Author SeaWulf settings. Both Python calls only reuse existing data.
#SBATCH --job-name=sinabro-final-checks
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/extra-final-checks-%j.out
#SBATCH --error=logs/extra-final-checks-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch from the project root}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
HUMAN="${SINABRO_EXTRA_HUMAN:-${PROJECT}/runs/aggregation_2305071}"
OUTPUT="${SINABRO_FINAL_CHECK_OUTPUT:-${PROJECT}/results/extra_analysis/baselines_stability_${SLURM_JOB_ID}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
CHECKS="${ANALYSIS_DIR}/Extra_analysis/03_baselines_and_stability"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); exact baselines and saved-path precision"
"$PYTHON_BIN" "$CHECKS/analyze_baselines.py" \
  --human-export "$HUMAN" \
  --signatures "$PACKAGE/sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt" \
  --out-dir "$OUTPUT/one_step"

"$PYTHON_BIN" "$CHECKS/analyze_stability.py" \
  --human-export "$HUMAN" --out-dir "$OUTPUT/stability" \
  --prefixes 100 200 500 1000 --per-bin 20 --sample-seed 42

echo "Both checks finished: $OUTPUT"
