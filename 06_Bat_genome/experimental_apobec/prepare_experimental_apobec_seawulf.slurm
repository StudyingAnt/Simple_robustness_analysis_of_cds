#!/usr/bin/env bash
# Source XLSX is transferred privately once; no experimental values are in Git.
#SBATCH --job-name=sinabro-apobec-prepare
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/apobec-prepare-%j.out
#SBATCH --error=logs/apobec-prepare-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
WORKBOOK="${SINABRO_APOBEC_WORKBOOK:-${PROJECT}/data/raw/bat_apobec/C to T Mutation data of six bat enzymes.xlsx}"
PROFILES="${SINABRO_APOBEC_PROFILES:-${PROJECT}/data/processed/experimental_apobec_v1}"
INPUTS="${SINABRO_APOBEC_INPUTS:-${PROJECT}/data/processed/experimental_apobec_inputs_v1}"
HUMAN="${SINABRO_HUMAN_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
BATS="${SINABRO_BAT_PROCESSED:-${PROJECT}/data/processed/bat1k_historical_v1}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); preparing experimental APOBEC inputs"
if [[ ! -f "$WORKBOOK" ]]; then
  echo "Missing original workbook: $WORKBOOK" >&2
  exit 1
fi
if [[ ! -e "$PROFILES" ]]; then
  "$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/prepare_profiles.py" \
    --workbook "$WORKBOOK" --out-dir "$PROFILES"
fi
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/prepare_inputs.py" \
  --human-input "$HUMAN" --bat-root "$BATS" --profiles-dir "$PROFILES" --out-dir "$INPUTS"
