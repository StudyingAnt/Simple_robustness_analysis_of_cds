# Manuscript and supplementary execution: author HPC record

The authors ran this stage on Stony Brook University's **SeaWulf** HPC system, which requires institutional access. For execution on your own computer, start with the [stage README](../README.md). These are the authors' execution examples, not a requirement for readers. If you use your own HPC system, adapt the scheduler, site/environment settings and paths before use.

PCA/features/statistics were computed in the combined downstream job; there was no separate stage-04 Slurm job. These scripts and reports describe execution on SeaWulf, rather than providing access to it.

## Execution guides and records

- [Actual downstream statistics and signature-sampler submissions](MAIN_ANALYSIS.md)
- [Matched archived S11 estimates; job 2345768](S11_HISTORICAL.md)
- [Fresh BRCA1 graph; job 2346048](BRCA1_GRAPH.md)
- [Actual observed-SNV annotation report and screenshot review](BRCA1_OBSERVED_SNVS.md)
- [Detailed recorded validation comparisons](records/DEVELOPMENT_VALIDATION_RECORD.md)
- [Historical supplementary source recovery](records/SUPPLEMENTARY_SOURCE_AUDIT.md)

## Retained Slurm files

- [analyze_s11_historical_seawulf.slurm](analyze_s11_historical_seawulf.slurm)
- [annotate_brca1_graph_seawulf.slurm](annotate_brca1_graph_seawulf.slurm)
- [run_brca1_graph_seawulf.slurm](run_brca1_graph_seawulf.slurm)
- [analyze_robustness_results_seawulf.slurm](analyze_robustness_results_seawulf.slurm)
- [validate_signature_sampler_seawulf.slurm](validate_signature_sampler_seawulf.slurm)

Commands in the execution records are issued from the authors' project root. Full raw/checkpoint data are external to Git. See the [author archive index](../../author_hpc/README.md) for the site/environment context and layout change.
