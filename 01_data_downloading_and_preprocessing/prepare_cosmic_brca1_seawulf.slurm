#!/usr/bin/env bash
# Verify the downloaded archive and inspect every BRCA1 source row; no login/network.
#SBATCH --job-name=cosmic-brca1-input
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/cosmic-brca1-input-%j.out
#SBATCH --error=logs/cosmic-brca1-input-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
ARCHIVE="${COSMIC_MUTANTCENSUS_TAR:-${PROJECT}/data/raw/cosmic/v101/GRCh38/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar}"
OUTPUT="${COSMIC_BRCA1_OUTPUT:-${PROJECT}/data/processed/cosmic_v101_grch38/brca1_source_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); verifying COSMIC archive and extracting BRCA1 source rows"
"$PYTHON_BIN" "$ANALYSIS_DIR/01_data_downloading_and_preprocessing/prepare_cosmic_brca1.py" \
  --archive "$ARCHIVE" --out-dir "$OUTPUT"
