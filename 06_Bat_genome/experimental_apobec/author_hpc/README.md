# Experimental APOBEC execution: author HPC record

The authors ran this stage on Stony Brook University's **SeaWulf** HPC system, which requires institutional access. For execution on your own computer, start with the [stage README](../README.md). These are the authors' execution examples, not a requirement for readers. If you use your own HPC system, adapt the scheduler, site/environment settings and paths before use.

Preparation job 2329085, serial/parallel pilot 2329255, complete audit 2336448 and export 2338741. The full campaign contained 765,672,000 paths. The author allocation controller remains at its original source path one level above; the direct reader runner does not require that controller or its pilot report.

## Execution guides and records

- [Actual preparation, pilot, allocation, audit and export procedure](EXECUTION.md)

## Retained Slurm files

- [export_experimental_results_seawulf.slurm](export_experimental_results_seawulf.slurm)
- [prepare_experimental_apobec_seawulf.slurm](prepare_experimental_apobec_seawulf.slurm)
- [run_experimental_apobec_seawulf.slurm](run_experimental_apobec_seawulf.slurm)
- [validate_experimental_pilot_seawulf.slurm](validate_experimental_pilot_seawulf.slurm)
- [validate_experimental_results_seawulf.slurm](validate_experimental_results_seawulf.slurm)

Commands in the execution records are issued from the authors' project root. Full raw/checkpoint data are external to Git. See the [author archive index](../../../author_hpc/README.md) for the site/environment context and layout change.
