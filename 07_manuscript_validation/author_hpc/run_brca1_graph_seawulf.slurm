#!/usr/bin/env bash
# Fresh BRCA1 graph: 1000 SBS2/BLOSUM80 paths. Observed COSMIC annotation is a later stage.
#SBATCH --job-name=sinabro-brca1-graph
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/brca1-graph-%j.out
#SBATCH --error=logs/brca1-graph-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPUs}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
PACKAGE="${SINABRO_SOURCE:-${PROJECT}/sinabro}"
INPUT="${SINABRO_FASTA:-${PROJECT}/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
OUTPUT="${SINABRO_BRCA1_OUTPUT:-${PROJECT}/runs/brca1_sbs2_blosum20_n1000_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); BRCA1 SBS2/BLOSUM80 graph, 1000 fresh paths"
"$PYTHON_BIN" "$ANALYSIS_DIR/07_manuscript_validation/run_supplementary_cases.py" \
  --stage s12 --input "$INPUT" --sinabro-source "$PACKAGE" --out-dir "$OUTPUT" \
  --n-traj 1000 --workers "$SLURM_CPUS_PER_TASK" --seed 42 --max-steps 10000
