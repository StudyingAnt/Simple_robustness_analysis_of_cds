#!/usr/bin/env bash
# All three datasets, all six profiles; serial versus eight workers.
#SBATCH --job-name=sinabro-apobec-pilot
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/apobec-pilot-%j.out
#SBATCH --error=logs/apobec-pilot-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
INPUTS="${SINABRO_APOBEC_INPUTS:-${PROJECT}/data/processed/experimental_apobec_inputs_v1}"
OUTPUT="${SINABRO_APOBEC_PILOT_OUTPUT:-${PROJECT}/runs/experimental_apobec_pilot_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); experimental APOBEC serial/parallel pilot"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/validate_experimental_pilot.py" \
  --inputs-dir "$INPUTS" --out-dir "$OUTPUT" --sinabro-source "$PACKAGE" \
  --workers "$SLURM_CPUS_PER_TASK" --n-traj 100 --seed 42 --max-steps 10000
