# Export audited results for manuscript analysis

Use aggregate_results.py after validate_results.py has passed. This stage exports the stored, independently checked summary statistics. It does not simulate trajectories, fit PCA, draw figures or test manuscript conclusions. NumPy and Biopython from the existing pinned environment are sufficient; pandas and a Sinabro import are not required.

## Audit and export on your computer

Run from the analysis repository root after the [full human calculation](README.md#reproduce-the-manuscript-calculation) finishes. Set `PROJECT` to your own absolute data/output directory. No scheduler or pilot report is required.

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
WORKERS=1

python 02_run_sinabro/validate_results.py \
  --run-dir "$PROJECT/runs/human_cosmic_n1000_seed42" \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --report-dir "$PROJECT/runs/human_audit" \
  --workers "$WORKERS" --expected-records 73214

python 02_run_sinabro/aggregate_results.py \
  --run-dir "$PROJECT/runs/human_cosmic_n1000_seed42" \
  --validation-dir "$PROJECT/runs/human_audit" \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --out-dir "$PROJECT/runs/human_export" \
  --signature-groups 02_run_sinabro/manuscript_signature_groups.json \
  --workers "$WORKERS" --expected-records 73214
```

Inspect `human_audit/validation_summary.json` for `validation_passed=true` and `all_paths_complete=true`, then `human_export/aggregation_manifest.json` for `aggregation_passed=true`. These directories are under `PROJECT/runs/`. Continue with [distribution figures](../03_mean_std_analysis/README.md), using `human_export` as the aggregation input.

The validation directory must contain a passing validation_summary.json and transcript_audit.csv from the same simulation. The manifest fingerprint, input hash, validated transcript set, dimensions and completed-path counts must agree. Every checkpoint's current SHA256 must match its audit-inventory entry before its summary arrays can be exported. Paths can change when transferring a run to another computer; matching content hashes establish correspondence. No unchecked or incomplete transcript is silently omitted.

The input must have GENCODE-style headers containing transcript ID, gene ID, transcript name and gene name. --input defaults to the path recorded in the audit report. --workers defaults to one and needs no scheduler variables. --expected-records is optional for other input sizes. For a small example using only a few signatures, omit --signature-groups; the manuscript subset file intentionally requires all its members to be present.

The output directory must be new and outside both the simulation and validation directories. Tables are built in a temporary sibling directory and published only after success. A mismatch stops export with a nonzero exit code and leaves no completed output directory. Errors appear in terminal output. Keep source results, inputs and audit files unchanged during export. Pick another output directory when repeating an export.

## Output tables and identifiers

All matrices preserve the simulation manifest's signature order. Rows are sorted by versioned transcript ID, with the same row correspondence across all tables. CSV files have no added row-index column.

| File | Meaning of rows and values |
| --- | --- |
| rums_mean_by_transcript_id.csv | Transcript column is the versioned ENST identifier; values are mean accepted synonymous steps per trajectory (RUMS). |
| rums_sd_by_transcript_id.csv | Same ENST rows; sample SD of trajectory counts within each transcript/signature (ddof=1). |
| rums_median_by_transcript_id.csv | Same ENST rows; median trajectory count. |
| rums_mc_se_by_transcript_id.csv | Same ENST rows; Monte Carlo SE of each mean, path SD / sqrt(number of paths). |
| all_gencode_rums_profile.csv | Same means, using the historical Transcript column: GENCODE transcript name, such as SPTBN1-201. This is the format expected by the existing plotting and PCA scripts. |
| all_gencode_rums_profile_noseqerr.csv | Historical-name means for the preserved 60-signature manuscript subset, when requested. |
| all_gencode_rums_profile_aetiology.csv | Historical-name means for the preserved 41-signature manuscript subset, when requested. |
| transcript_annotations.csv | Explicit ENST/name/gene-ID/gene-name mapping, CDS length, internal-stop flag and complete FASTA identifier. |
| signature_summary.csv | Descriptive mean, sample SD, median, min and max of transcript-level RUMS, plus median Monte Carlo SE and descending ranks. |
| aggregation_manifest.json | Source/configuration and audit hashes, output hashes, row/signature counts, grouping/weighting rules, versions and script hash. |

The full manuscript matrices contain 73,214 rows and 79 signature columns, plus the Transcript column. The two optional subset matrices have 60 and 41 signature columns. transcript_annotations.csv preserves all 73,214 records; there is no representative-transcript selection or gene-level averaging. Local inspection of the verified GENCODE input found 19,309 distinct versioned gene IDs and 73 internal-stop transcripts; both totals are reported again during export. Internal-stop records remain included, matching the simulation/preprocessing policy.

The historical simulator named its CSV files using FASTA field 5 (transcript name), whereas run_rums.py names NPZ files using field 1 (ENST). Both representations are exported deliberately, with a one-to-one mapping. Duplicate ENST IDs or duplicate transcript names cause an error rather than an ambiguous join. Names and versions are preserved verbatim.

manuscript_signature_groups.json copies membership and order from the preserved reduce_signatures.py; it does not reclassify signatures using a newer COSMIC release. Tests also compare these lists with run_pca.py. The all-signature column order always comes from the actual run manifest.

## Different sources of variation

- rums_sd_by_transcript_id.csv measures trajectory variability for one fixed transcript/signature pair. It is not the across-transcript SD used to rank signatures in the historical violin analysis.
- signature_summary.csv uses mean_rums_across_transcripts and sd_rums_across_transcripts, giving every transcript isoform equal weight. The sample SD uses ddof=1, matching R's sd() for the same values. Transcript isoforms sharing a gene are retained; this is a descriptive summary, not a claim of biological independence.
- Monte Carlo SE measures simulation uncertainty in one transcript's estimated RUMS. It does not treat 1,000 trajectories as 1,000 biological replicates. With one trajectory, path SD and Monte Carlo SE are undefined (NaN); with one transcript, across-transcript SD/rank is undefined.
- Descending ranks use the signature name alphabetically to break exact ties. The summary includes all selected signatures; no significance tests or new filtering policies are applied.

## Connection to the existing manuscript workflow

| Existing script | Input supplied by this stage | Work remaining before executing that script |
| --- | --- | --- |
| 03_mean_std_analysis/analyze_signature_distributions.py | Completed export directory with canonical/legacy matrices, annotations, summary and manifest. | Use the portable [descriptive analysis workflow](../03_mean_std_analysis/README.md) to check ranks and render distributions in the pinned environment. |
| 03_mean_std_analysis/draw_violin_plots_by_rank.R (historical reference) | all_gencode_rums_profile.csv | Original paths/layout are preserved; use the Python workflow above for the current analysis. |
| 04_PCA_and_gc_contents/reduce_signatures.py | Its two outputs are already exported when the manuscript groups JSON is provided. | No need to rerun the historical reduction script. |
| 04_PCA_and_gc_contents/run_pca.py | The full, noseqerr and aetiology matrices, with the expected names and column order. | Make paths portable, record the PCA fitting/standardization policy and check outputs before plotting. |
| 04_PCA_and_gc_contents/compute_gc_contents_sequences.py and draw_pca_plots.R | transcript_annotations.csv provides the explicit ID/name correspondence. | Compute/check GC features from the same CDS-plus file and join by an explicit key; do not assume that independently sorted rows match. |

The historical combine_rums_profile.py is retained for its original trajectory-CSV workflow. Use aggregate_results.py for audited NPZ outputs. Numeric equality with old manuscript results is not assumed because these outputs come from the corrected Sinabro implementation. Keep this export under PROJECT/runs until checked; final analysis tables and figures belong under PROJECT/results.

## Author execution

The authors used SeaWulf for the full export. Their site-specific wrapper and original paths are retained under [author_hpc/](author_hpc/README.md); access to that system is not required to use the portable exporter.
