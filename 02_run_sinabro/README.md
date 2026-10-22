# Human RUMS calculation

Use the [installed environment and Sinabro package](../CONDA_SETUP.md). All commands below run from the analysis repository root in a terminal; no scheduler or institutional account is required. Set `PROJECT` to an absolute directory for your inputs and generated outputs, and `WORKERS` to the CPU processes you want to use.

## Try the workflow on the included example

The [example FASTA](example_cdsplus.fa) contains three **synthetic** CDS-plus records with GENCODE-style identifiers. It is an execution example, not a human reference or a manuscript result. Choose a new output directory.

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
WORKERS=1

python 02_run_sinabro/run_rums.py \
  --input 02_run_sinabro/example_cdsplus.fa \
  --out-dir "$PROJECT/runs/example_rums" \
  --signatures SBS1 SBS2 --n-traj 20 --seed 42 --workers "$WORKERS"

python 02_run_sinabro/validate_results.py \
  --run-dir "$PROJECT/runs/example_rums" \
  --input 02_run_sinabro/example_cdsplus.fa \
  --report-dir "$PROJECT/runs/example_audit" --workers "$WORKERS"

python 02_run_sinabro/aggregate_results.py \
  --run-dir "$PROJECT/runs/example_rums" \
  --validation-dir "$PROJECT/runs/example_audit" \
  --input 02_run_sinabro/example_cdsplus.fa \
  --out-dir "$PROJECT/runs/example_export" --workers "$WORKERS"
```

This calculates 120 paths (3 records × 2 signatures × 20). The audit should report `validation_passed=true` and `all_paths_complete=true`; export should report `aggregation_passed=true`. Continue with [distribution plots](../03_mean_std_analysis/README.md) using `example_export` as the input. The example does not require downloading GENCODE.

## Reproduce the manuscript calculation

First [download GENCODE v40 and prepare CDS-plus](../01_data_downloading_and_preprocessing/README.md). With the same `PROJECT` and `WORKERS` variables:

```bash
python 02_run_sinabro/run_rums.py \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --out-dir "$PROJECT/runs/human_cosmic_n1000_seed42" \
  --n-traj 1000 --seed 42 --max-steps 10000 \
  --workers "$WORKERS" --expected-records 73214
```

Omitting `--signatures` selects all 79 COSMIC v3.3.1 GRCh37 SBS profiles. The full calculation contains **5,783,906,000 paths** and can take substantial time on a personal computer. The commands specify the calculation; they do not imply that it is a short workstation run. Smaller trajectory counts or input subsets are examples rather than full manuscript reproduction.

The stopping criterion is the first amino-acid sequence change. RUMS counts accepted neutral steps before that event. The 10,000-step cap is a safety limit; capped or unavailable-event paths are incomplete and must not be interpreted as complete robustness estimates.

Audit and export the completed full calculation using [AGGREGATION.md](AGGREGATION.md). Then follow [distribution analysis](../03_mean_std_analysis/README.md) and [manuscript analyses](../07_manuscript_validation/README.md). [REPRODUCTION.md](REPRODUCTION.md) explains input format, resume and execution settings.

## How the authors ran this stage

The authors performed the full calculation on Stony Brook University's SeaWulf HPC system. SeaWulf requires institutional access and is not required for these Python commands. The actual Slurm wrappers, run records and serial/parallel pilot materials are retained under [author_hpc/](author_hpc/README.md). Researchers using their own HPC system may adapt those examples to their site's scheduler, partitions, environment and paths.
