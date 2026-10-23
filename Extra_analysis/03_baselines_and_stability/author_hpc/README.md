# Author SeaWulf execution

These are the authors' institutional paths and scheduler settings. Adapt them
for another system or use the two portable Python commands in the parent README.

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/Extra_analysis/03_baselines_and_stability/author_hpc/run_baselines_stability_seawulf.slurm
```

One job sequentially runs both checks, with one CPU and 8 GB. The three-hour
limit is a resource request, not a measured runtime estimate. The full raw human
run and audit recorded in `runs/aggregation_2305071/aggregation_manifest.json`
must remain available for the second check. No mutation simulation is run.

Defaults may be changed with `PROJECT`, `ANALYSIS_DIR`, `SINABRO_ENV`,
`SINABRO_SOURCE`, `SINABRO_EXTRA_HUMAN` and `SINABRO_FINAL_CHECK_OUTPUT`.
For relocated raw runs/audits, use the Python `--run-dir`/`--audit-dir` options;
the input checksums must remain identical.

After the job, substitute its job number and inspect both reports:

```bash
OUT=results/extra_analysis/baselines_stability_JOBID
cat "$OUT/one_step/analysis_summary.json"
cat "$OUT/stability/analysis_summary.json"
```

Both should have `analysis_passed: true` and `fresh_simulation: false`. Baseline
production output should have `bounded_subset: false`, 19,309 representatives and
79 signatures. Stability should have `bounded_export: false`, prefixes
100/200/500/1000, at most 184 metadata-selected representatives and all 79 profiles.
It deliberately records `stability_sample_is_full_genome: false` and does not
declare universal sufficiency of 1,000 paths.
