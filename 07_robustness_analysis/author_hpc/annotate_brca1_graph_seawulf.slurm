#!/usr/bin/env bash
# Annotate existing audited graph/data; no new trajectories or downloads.
#SBATCH --job-name=brca1-snv-matches
#SBATCH --partition=short-96core-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/brca1-snv-matches-%j.out
#SBATCH --error=logs/brca1-snv-matches-%j.err
#SBATCH --mail-user=joon-hyun.song@stonybrook.edu
#SBATCH --mail-type=END,FAIL
set -eo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch}"
PROJECT="${PROJECT:-${SLURM_SUBMIT_DIR}}"
ANALYSIS_DIR="${ANALYSIS_DIR:-${PROJECT}/Simple_robustness_analysis_of_cds}"
SINABRO_ENV="${SINABRO_ENV:-${PROJECT}/sinabro-revision}"
OUTPUT="${COSMIC_GRAPH_OUTPUT:-${PROJECT}/results/brca1_cosmic_snv_matches_${SLURM_JOB_ID}}"
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SINABRO_ENV"
set -u
PYTHON_BIN="${SINABRO_PYTHON:-$CONDA_PREFIX/bin/python}"
if [[ -n "${SINABRO_BRCA1_GRAPH_DIR:-}" ]]; then
  GRAPH_ARGS=(--graph-dir "$SINABRO_BRCA1_GRAPH_DIR")
else
  GRAPH_ARGS=(--graph-root "$PROJECT/runs")
fi
if [[ -n "${COSMIC_BREAST_DIR:-}" ]]; then
  VARIANT_ARGS=(--variants-dir "$COSMIC_BREAST_DIR")
else
  VARIANT_ARGS=(--variants-root "$PROJECT/data/processed/cosmic_v101_grch38")
fi
cd -- "$PROJECT"
echo "Job $SLURM_JOB_ID on $(hostname); matching COSMIC breast SNVs to the existing BRCA1 graph"
"$PYTHON_BIN" "$ANALYSIS_DIR/07_robustness_analysis/annotate_brca1_graph.py" \
  "${GRAPH_ARGS[@]}" "${VARIANT_ARGS[@]}" --out-dir "$OUTPUT"
