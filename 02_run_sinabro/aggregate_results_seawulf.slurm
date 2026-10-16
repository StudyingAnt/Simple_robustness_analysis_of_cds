#!/usr/bin/env bash
# Export a previously validated run; submit from the HPC project root.
#SBATCH --job-name=sinabro-aggregate
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/aggregate-%j.out
#SBATCH --error=logs/aggregate-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL

set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
: "${SINABRO_AGGREGATE_VALIDATION:?Export the passing results_check_JOBID directory before sbatch}"

PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
INPUT_FASTA="${SINABRO_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
RUN_DIR="${SINABRO_AGGREGATE_RUN:-${PROJECT}/runs/hpc_gencode_v40_n1000_seed42}"
OUTPUT_DIR="${SINABRO_AGGREGATE_OUTPUT:-${PROJECT}/runs/aggregation_${SLURM_JOB_ID}}"

module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); exporting audited results: $RUN_DIR"
"$PYTHON_BIN" "$ANALYSIS_DIR/02_run_sinabro/aggregate_results.py" \
  --run-dir "$RUN_DIR" \
  --validation-dir "$SINABRO_AGGREGATE_VALIDATION" \
  --input "$INPUT_FASTA" \
  --out-dir "$OUTPUT_DIR" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --signature-groups "$ANALYSIS_DIR/02_run_sinabro/manuscript_signature_groups.json" \
  --expected-records 73214
