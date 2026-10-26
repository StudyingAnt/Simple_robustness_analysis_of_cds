#!/usr/bin/env bash
# Run from the analysis repository root in the pinned Conda environment.
set -euo pipefail
: "${HUMAN_EXPORT:?Set the authenticated aggregation directory}"
: "${HUMAN_FASTA:?Set the relocated original CDS-plus FASTA}"
: "${SINABRO_SOURCE:?Set the Sinabro 3.0.0 checkout directory}"
: "${OUT_ROOT:?Set a new output root}"
PYTHON="${PYTHON:-python}"
WORKERS="${WORKERS:-8}"
ANALYSIS=Extra_analysis/04_composition_and_codon_null
mkdir -p "$OUT_ROOT"
"$PYTHON" "$ANALYSIS/prepare_inputs.py" \
  --human-export "$HUMAN_EXPORT" --human-fasta "$HUMAN_FASTA" \
  --signatures "$SINABRO_SOURCE/sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt" \
  --out-dir "$OUT_ROOT/inputs"
"$PYTHON" "$ANALYSIS/analyze_composition.py" \
  --inputs "$OUT_ROOT/inputs" --out-dir "$OUT_ROOT/composition"
"$PYTHON" "$ANALYSIS/run_codon_null.py" \
  --inputs "$OUT_ROOT/inputs" --sinabro-source "$SINABRO_SOURCE" \
  --out-dir "$OUT_ROOT/null_n500_r20" --n-traj 500 --recodings 20 --workers "$WORKERS"
"$PYTHON" "$ANALYSIS/analyze_codon_null.py" \
  --inputs "$OUT_ROOT/inputs" --simulation "$OUT_ROOT/null_n500_r20" \
  --out-dir "$OUT_ROOT/codon_null"
