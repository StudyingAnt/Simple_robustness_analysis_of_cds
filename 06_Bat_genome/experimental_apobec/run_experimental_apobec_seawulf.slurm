#!/usr/bin/env bash
# One complete submission, 50 allocation tasks, 100 scientific shards/dataset.
#SBATCH --job-name=sinabro-apobec
#SBATCH --partition=short-96core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --array=0-49%8
#SBATCH --time=04:00:00
#SBATCH --output=logs/apobec-%A_%a.out
#SBATCH --error=logs/apobec-%A_%a.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_ARRAY_TASK_ID:?Requires an array}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
: "${SINABRO_APOBEC_PILOT:?Set the passing experimental_apobec_pilot_JOBID directory}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
INPUTS="${SINABRO_APOBEC_INPUTS:-${PROJECT}/data/processed/experimental_apobec_inputs_v1}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
SEED="${SINABRO_SEED:-42}"
OUTPUT="${SINABRO_APOBEC_RUN:-${PROJECT}/runs/hpc_experimental_apobec_n${N_TRAJ}_seed${SEED}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID on $(hostname); complete experimental APOBEC campaign"
"$PYTHON_BIN" "$ANALYSIS_DIR/06_Bat_genome/experimental_apobec/run_campaign.py" \
  --inputs-dir "$INPUTS" --pilot-dir "$SINABRO_APOBEC_PILOT" --out-root "$OUTPUT" \
  --task-index "$SLURM_ARRAY_TASK_ID" --task-count 50 --shards-per-dataset 100 \
  --sinabro-source "$PACKAGE" --workers "$SLURM_CPUS_PER_TASK" --n-traj "$N_TRAJ" --seed "$SEED" --max-steps 10000
