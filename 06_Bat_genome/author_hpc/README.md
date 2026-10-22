# Bat extraction and COSMIC execution: author HPC record

The authors ran this stage on Stony Brook University's **SeaWulf** HPC system, which requires institutional access. For execution on your own computer, start with the [stage README](../README.md). These are the authors' execution examples, not a requirement for readers. If you use your own HPC system, adapt the scheduler, site/environment settings and paths before use.

Author pilot job 2324062; complete audit job 2327882 and export job 2328086. The frozen populations are 31,185 and 23,213 transcripts; the full campaign contained 435,184,000 paths. The allocation controller retains its original source path one level above for reproducible provenance.

## Execution guides and records

- [Raw source download, CDS extraction and author pilot](PREPARATION.md)
- [Full array, audit and export](COSMIC_EXECUTION.md)
- [Recorded cohorts and completed-run overview](REFERENCE_RUN.md)
- [Local extraction and bounded validation record](records/BAT_REFERENCE_CHECKS.md)

## Retained Slurm files

- [aggregate_bat_results_seawulf.slurm](aggregate_bat_results_seawulf.slurm)
- [download_bat_inputs_seawulf.slurm](download_bat_inputs_seawulf.slurm)
- [extract_bat_cds_seawulf.slurm](extract_bat_cds_seawulf.slurm)
- [run_bat_cosmic_seawulf.slurm](run_bat_cosmic_seawulf.slurm)
- [validate_bat_pilot_seawulf.slurm](validate_bat_pilot_seawulf.slurm)
- [validate_bat_results_seawulf.slurm](validate_bat_results_seawulf.slurm)

Commands in the execution records are issued from the authors' project root. Full raw/checkpoint data are external to Git. See the [author archive index](../../author_hpc/README.md) for the site/environment context and layout change.
