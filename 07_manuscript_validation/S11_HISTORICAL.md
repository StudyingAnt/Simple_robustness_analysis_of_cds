# S11: matched-transcript comparison of historical estimates

This stage uses bundled archived estimates for the same 3,740 transcripts under protein identity and BLOSUM80 (threshold 20). It does not simulate new paths. The pinned human FASTA verifies transcript identity; original HDF5 files are not needed.

From the analysis repository root with the environment active:

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
python 07_manuscript_validation/analyze_s11_historical.py \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --out-dir "$PROJECT/results/s11_historical"
```

Choose a new output directory. Inspect `s11_historical_summary.json`: expected fields include `analysis_passed=true`, `historical_results_only=true`, `n_transcripts=3740` and `input_sequences_verified=3740`. Outputs include matched transcript values, criterion distributions and a scatter plot (CSV/PDF/PNG).

Reference means of the archived estimates are **0.8116874332** (identity) and **3.8720010695** (BLOSUM). The original plotted populations were 73,214 and 3,740; this workflow compares the same 3,740 transcripts. Historical seeds and per-path counts were not recovered. Matched transcripts do not imply paired random paths, a Monte Carlo confidence interval or a fresh corrected-package BLOSUM campaign. BLOSUM distance is an alignment proxy rather than measured functional loss.

The authors ran this postprocessing on SeaWulf. Their [submitted script and complete execution record](author_hpc/S11_HISTORICAL.md) are preserved separately; institutional access is not required for the command above. Researchers using their own HPC may adapt the example to their site's settings. [Source recovery](author_hpc/records/SUPPLEMENTARY_SOURCE_AUDIT.md) documents archived provenance.
