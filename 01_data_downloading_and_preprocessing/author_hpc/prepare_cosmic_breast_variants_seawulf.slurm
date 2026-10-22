#!/usr/bin/env bash
# Join COSMIC phenotypes and map breast BRCA1 SNVs to the existing simulated CDS.
#SBATCH --job-name=cosmic-breast-map
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/cosmic-breast-map-%j.out
#SBATCH --error=logs/cosmic-breast-map-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
RAW_DIR="${COSMIC_RAW_DIR:-${PROJECT}/data/raw/cosmic/v101/GRCh38}"
INPUT="${SINABRO_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
OUTPUT="${COSMIC_BREAST_OUTPUT:-${PROJECT}/data/processed/cosmic_v101_grch38/brca1_breast_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); COSMIC breast selection and BRCA1 target-CDS mapping"
"$PYTHON_BIN" "$ANALYSIS_DIR/01_data_downloading_and_preprocessing/prepare_cosmic_breast_variants.py" \
  --raw-dir "$RAW_DIR" --input "$INPUT" --out-dir "$OUTPUT"
