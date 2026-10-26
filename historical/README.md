# Historical analysis sources

This directory preserves earlier analysis and plotting sources, exploratory examples. Current execution starts in the [manuscript reproduction guide](../MANUSCRIPT_REPRODUCTION.md). These archived scripts may use old package imports, fixed paths, SGE or a separate R environment; they are method provenance rather than the current command-line interface.

| Directory | Preserved material |
| --- | --- |
| `02_run_sinabro/` | Earlier single-gene runner, CSV combiner and SGE wrapper |
| `03_mean_std_analysis/`, `04_PCA_and_gc_contents/` | Earlier plotting, PCA and sequence-feature scripts |
| `05_Cancer_census_genes/` | Earlier Cancer Gene Census extraction sources |
| `06_Bat_genome/` | Earlier bat extraction, combination and R plotting scripts |
| `99_Supplementary/` | Earlier supporting analyses and graph plotting |
| `examples/` | Earlier bounded supplementary Slurm example |

The original analysis snapshot is commit `c2c271e36a0f478fd3b3464902a450c514f5b789`. The complete pre-publication checkout is preserved on [backup/before-paper-release-2026-10-21](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-paper-release-2026-10-21). The preceding public main is preserved separately on [backup/main-before-paper-release-2026-10-21](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/main-before-paper-release-2026-10-21).

The frozen signature-subset manifest retains its original `04_PCA_and_gc_contents/reduce_signatures.py` source path and checksum. That small source remains at its original location for independent membership/order checks. Current PCA is implemented in `07_robustness_analysis/analyze_robustness_results.py`.

Author pilot, source-recovery and completed-run validation records are retained inside each stage’s `author_hpc/records/` folder; see the [author execution index](../author_hpc/README.md).
