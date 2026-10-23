# Additional analyses

New analyses beyond the original-paper reproduction are kept here. The numbered
manuscript workflows, original outputs and Sinabro package are not modified.
Future additional/reviewer analyses should use another clearly named subfolder
here and write results to a new external output directory.

| Analysis | Inputs | Main outputs |
|---|---|---|
| [01_statistics_and_sequence_features](01_statistics_and_sequence_features/README.md) | Audited exported means, authenticated annotations and the original CDS-plus FASTAs | Human Friedman/paired Wilcoxon results, paired effect sizes and median CIs, three-cohort GC3 distributions, outlier annotations, descriptive species comparisons |
| [02_bat_assay_figure](02_bat_assay_figure/README.md) | Existing authenticated six-enzyme motif means | Figure 6 with a common percentage axis, draft caption and provenance; no new assay or RUMS calculation |

Start with the portable instructions in the analysis folder. The authors' SeaWulf
submission example is stored separately in `author_hpc/`; it requires the authors'
institutional access and paths. Other users can run Python directly or adapt the
resource/environment settings for their own scheduler.

No new RUMS trajectories are needed for the first analysis. A passing execution
means that the stated postprocessing and input checks completed; it does not
establish adaptation, coevolution, functional loss or a lifespan association.
