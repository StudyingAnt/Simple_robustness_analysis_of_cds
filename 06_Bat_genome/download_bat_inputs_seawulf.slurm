#!/usr/bin/env bash
# Raw reference download only; independent of human downstream analyses.
#SBATCH --job-name=sinabro-bat-download
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/bat-download-%j.out
#SBATCH --error=logs/bat-download-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
OUTPUT_DIR="${SINABRO_BAT_RAW:-${PROJECT}/data/raw/bat1k_pilot_2020}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); downloading raw Bat1K inputs: $OUTPUT_DIR"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/download_bat_inputs.py" --out-dir "$OUTPUT_DIR"
