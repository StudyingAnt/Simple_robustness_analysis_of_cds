#!/usr/bin/env bash
# Explicit supplementary reconstructions; source COSMIC variant overlay is unavailable.
#SBATCH --job-name=sinabro-supplement
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/paper-supplement-%j.out
#SBATCH --error=logs/paper-supplement-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
STAGE="${SINABRO_SUPPLEMENT:-s11}"
case "$STAGE" in s11) N_TRAJ=1000;; s12) N_TRAJ=100;; *) echo 'Use s11 or s12' >&2; exit 1;; esac
INPUT="${SINABRO_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
OUTPUT="${SINABRO_SUPPLEMENT_OUTPUT:-${PROJECT}/results/supplementary_${STAGE}_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); supplementary $STAGE reconstruction"
"$PYTHON_BIN" "$ANALYSIS_DIR/07_robustness_analysis/run_supplementary_cases.py" \
  --stage "$STAGE" --input "$INPUT" --sinabro-source "$PACKAGE" --out-dir "$OUTPUT" \
  --n-traj "$N_TRAJ" --workers "$SLURM_CPUS_PER_TASK" --seed 42 --max-steps 10000
