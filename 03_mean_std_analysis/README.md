# Descriptive RUMS ranks and distribution figures

Run this stage after the [audited table export](../02_run_sinabro/AGGREGATION.md) completes. analyze_signature_distributions.py reads the completed export, checks its integrity and descriptive statistics, then draws distribution figures. It uses the existing Python/NumPy/SciPy/Matplotlib pinned environment; no R packages, Sinabro import, simulation or new environment are needed.

## Portable execution

From the analysis repository root:

```bash
python 03_mean_std_analysis/analyze_signature_distributions.py \
  --aggregation-dir /PATH/TO/PROJECT/runs/human_export \
  --out-dir /PATH/TO/PROJECT/results/signature_distributions \
  --workers 1 --expected-records 73214 --expected-signatures 79
```

The output directory must be new and outside the aggregation directory. Intermediate files are staged separately and the completed directory is published only after success. A failure returns a nonzero exit code and leaves the original export unchanged. Keep the input export fixed during execution.

For a smaller audited export, omit the two expected-dimension assertions. Defaults: one worker, 15 signatures per page, 256 density evaluation points and PNG resolution 200 dpi. --workers, --page-size, --density-points (at least 32) and --dpi can be supplied explicitly. The Python code has no scheduler paths or required Slurm variables. The saved aggregation manifest supplies the signature membership/order and dataset dimensions.

## Check the output

Read `analysis_summary.json` in your chosen output directory. It should report `analysis_passed=true` and `summary_and_rank_agreement=true`. Full human data should report 73,214 transcripts and 79 signatures. For the small included example, use its `example_export` directory and omit `--expected-records` and `--expected-signatures`; two signatures and three synthetic transcripts are not manuscript results.

## Author execution

The authors ran the full-data analysis on Stony Brook University's SeaWulf system, which requires institutional access. The [actual submission settings and run record](author_hpc/README.md) are retained separately. The portable command above needs no HPC account; adapt the optional Slurm example to your own site if needed.

## Checks and statistical definitions

The program checks the aggregation manifest/configuration dimensions and SHA256 of every listed input CSV. It verifies exact ENST/name correspondence using transcript_annotations.csv and compares historical-name mean matrices (including requested subsets) with the canonical ENST matrix. Gene counts and retained internal-stop counts must match the export manifest.

It recalculates each signature's across-transcript mean, sample SD (ddof=1), median, minimum and maximum from the canonical RUMS matrix. These must match signature_summary.csv to rtol=1e-10 and atol=1e-12, and the full-set mean/SD ranks must agree. The per-signature sample is the set of transcript RUMS estimates, with all isoforms equally weighted. SD here describes variability of mean RUMS between transcripts; it is distinct from the trajectory-count SD and Monte Carlo SE within one transcript/signature pair.

The all, noseqerr and aetiology groups use the 79/60/41 memberships preserved in the input export. Ranking within each group is descending, with exact ties resolved alphabetically by signature name. No gene-level collapse, internal-stop exclusion, transcript subsampling, p-values or independence assumptions are introduced. This stage does not rerun the historical Kruskal-Wallis/Dunn/Levene tests or establish a revised inferential result.

## Files and rendering

- signature_rankings.csv: full descriptive values and within-group mean/SD ranks, retaining every signature.
- all_rank_by_mean.pdf, all_rank_by_sd.pdf and corresponding group files: multipage PDFs of horizontal violins, with 15 signatures per page by default.
- GROUP_METRIC_START_END.png: PNG of every PDF rank page, using the same underlying data.
- aetiology_top_bottom.pdf/png: four-panel top/bottom comparison under mean and SD ranks; the full dataset uses 15 signatures per panel from the preserved 41-signature group.
- analysis_summary.json: observed leaders and SBS2 ranks, agreement checks, source/output hashes, environment and rendering settings.

Every violin uses all transcript values for that signature with a Gaussian KDE and Scott bandwidth. KDEs are calculated once per signature and reused across panels. Widths are normalized separately for each violin; median and interquartile range are overlaid. Constant-value columns are represented by a line instead of fitting a singular KDE. The shared RUMS axis spans the entire observed range, so values above the historical 6.5 limit remain visible. SBS2 is orange and other signatures blue; highlighting does not affect ranks or data selection.

These figures reproduce the descriptive quantities and ranked distribution comparisons, with new Python rendering. They do not reproduce the historical images pixel for pixel: ggplot2's KDE bandwidth and rendering differ, dense random jitter is replaced by median/IQR marks, fonts/layout are explicit, and the full range is retained. The original draw_violin_plots_by_rank.R remains preserved for reference. Its aetiology section overwrites its rank ranges and ultimately renders only bottom ranks 32-41; the Python workflow exports the complete ranking and top/bottom-15 comparison described in the manuscript.

PDFs are suitable for subsequent layout/editing; PNGs are previews. Validate full-run numerical results and inspect the exported figures before treating them as final manuscript panels. Plot file bytes can differ across runs because PDF metadata and library versions affect rendering; input hashes, numerical settings and output hashes are recorded for provenance.
