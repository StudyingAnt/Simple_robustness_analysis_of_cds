#!/usr/bin/env bash
# Submit from the HPC project root after the raw download job has passed.
#SBATCH --job-name=sinabro-bat-cds
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/bat-cds-%j.out
#SBATCH --error=logs/bat-cds-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
RAW_DIR="${SINABRO_BAT_RAW:-${PROJECT}/data/raw/bat1k_pilot_2020}"
PROCESSED_ROOT="${SINABRO_BAT_PROCESSED:-${PROJECT}/data/processed/bat1k_historical_v1}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); extracting frozen historical bat cohorts: $PROCESSED_ROOT"

for SPECIES in molMol2 myoMyo6; do
  "$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/extract_bat_cds.py" \
    --raw-dir "$RAW_DIR" \
    --out-dir "$PROCESSED_ROOT/$SPECIES" \
    --species "$SPECIES"
done
