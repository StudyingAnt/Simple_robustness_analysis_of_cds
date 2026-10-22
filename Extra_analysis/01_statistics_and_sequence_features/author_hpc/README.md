# Author SeaWulf execution

This submission script is for the authors' institutional SeaWulf account. Other
users should use the portable Python command in the parent README, or adapt the
partition, modules, environment, mail and storage settings for their own system.
This folder contains an execution recipe, not evidence that the full additional
analysis has already run.

The inputs below are the previously completed and audited exports. The job reads
them and writes a new results folder; it does not rerun RUMS.

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/Extra_analysis/01_statistics_and_sequence_features/author_hpc/run_statistics_seawulf.slurm
```

One job, one node, one CPU, 8 GB memory, three-hour time limit. This is a resource
limit, not a measured duration estimate. Defaults:

- Human: `runs/aggregation_2305071`
- Bat COSMIC: `runs/bat_cosmic_export_2328086`
- Experimental APOBEC: `runs/experimental_apobec_export_2338741`
- Environment: `sinabro-revision`
- Output: `results/extra_analysis/statistics_JOBID`

Paths are relative to the project root above. Override with `PROJECT`,
`ANALYSIS_DIR`, `SINABRO_ENV`, `SINABRO_EXTRA_HUMAN`, `SINABRO_EXTRA_BATS`,
`SINABRO_EXTRA_APOBEC` and `SINABRO_EXTRA_OUTPUT` if needed. The Python script can
also accept relocated, hash-identical FASTAs; see the parent README.

After completion, substitute the submitted job number:

```bash
cat results/extra_analysis/statistics_JOBID/additional_analysis_summary.json
```

Check `analysis_passed`, `bounded_subset: false`, 73,214 human transcripts,
19,309 representative genes, 79 signatures, 3,081 paired comparisons, and three
sequence-feature datasets. A changed annotation may alter the expected gene
count and should be investigated. `species_confidence_intervals_completed` is
deliberately false: verified bat gene grouping is still needed for that request.
