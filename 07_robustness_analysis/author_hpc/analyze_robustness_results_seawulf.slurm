#!/usr/bin/env bash
# Original-paper postprocessing only; existing simulated paths are reused.
#SBATCH --job-name=sinabro-paper-results
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=logs/paper-results-%j.out
#SBATCH --error=logs/paper-results-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
HUMAN="${SINABRO_PAPER_HUMAN:-${PROJECT}/runs/aggregation_2305071}"
BATS="${SINABRO_PAPER_BATS:-${PROJECT}/runs/bat_cosmic_export_2328086}"
APOBEC="${SINABRO_PAPER_APOBEC:-${PROJECT}/runs/experimental_apobec_export_2338741}"
OUTPUT="${SINABRO_PAPER_OUTPUT:-${PROJECT}/results/robustness_analysis_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); robustness comparisons, PCA and sequence-feature analysis"
"$PYTHON_BIN" "$ANALYSIS_DIR/07_robustness_analysis/analyze_robustness_results.py" \
  --human-export "$HUMAN" --bat-export "$BATS" --apobec-export "$APOBEC" \
  --signatures "$PACKAGE/sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt" \
  --out-dir "$OUTPUT" --n-traj 1000
