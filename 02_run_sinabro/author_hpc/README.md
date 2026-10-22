# Human RUMS execution: author HPC record

The authors ran this stage on Stony Brook University's **SeaWulf** HPC system, which requires institutional access. For execution on your own computer, start with the [stage README](../README.md). These are the authors' execution examples, not a requirement for readers. If you use your own HPC system, adapt the scheduler, site/environment settings and paths before use.

Reference full campaign: 73,214 transcripts × 79 signatures × 1,000 paths, 100 scientific shards. The full audit passed without incomplete paths; reference export job 2305071. Pilot files here are author verification material, not a required reader workflow.

## Execution guides and records

- [Full array, completion audit, resume and table export](EXECUTION.md)
- [Author serial/parallel execution checks](PARALLEL_CHECKS.md)
- [Recorded full-data exporter setup](EXPORT.md)
- [20-path/100-path timing and validation record](records/HUMAN_PILOTS.md)

## Retained Slurm files

- [aggregate_results_seawulf.slurm](aggregate_results_seawulf.slurm)
- [run_rums_seawulf.slurm](run_rums_seawulf.slurm)
- [validate_parallel_seawulf.slurm](validate_parallel_seawulf.slurm)
- [validate_results_seawulf.slurm](validate_results_seawulf.slurm)

Commands in the execution records are issued from the authors' project root. Full raw/checkpoint data are external to Git. See the [author archive index](../../author_hpc/README.md) for the site/environment context and layout change.
