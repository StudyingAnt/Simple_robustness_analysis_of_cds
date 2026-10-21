# Human RUMS simulation, audit and export

The manuscript configuration uses 73,214 GENCODE v40 CDS-plus transcripts, 79 COSMIC SBS signatures and 1,000 independent paths per transcript/signature combination. The stopping criterion is the first protein-sequence change. Completed paths report accepted neutral steps before that event; capped or otherwise incomplete paths cannot contribute a complete robustness value.

## Current workflow

1. Prepare the [pinned environment](../CONDA_SETUP.md) and [GENCODE input](../01_data_downloading_and_preprocessing/README.md).
2. Follow [portable execution and parallel validation](REPRODUCTION.md) to compare serial and parallel results on a small selection.
3. Use [HPC submission and resume](HPC_SUBMISSION.md) for the SeaWulf job array, or run `run_rums.py` directly on another computer.
4. Audit every saved archive with `validate_results.py` before [exporting tables](AGGREGATION.md).
5. Use [signature distributions](../03_mean_std_analysis/README.md) and [manuscript downstream analyses](../07_manuscript_validation/README.md) to generate numerical tables and plots from the checked exports.

Keep raw references in project `data/raw/`, sequence inputs in `data/processed/`, parameterized simulations/checkpoints and intermediate exports in `runs/`, and final analysis tables and plots in `results/`. All paths are configurable. Python runners are portable; `_seawulf.slurm` files contain cluster-specific settings.

The completed human audit and export covered all 73,214 × 79 × 1,000 paths without incomplete paths, according to the supplied HPC logs. This establishes checkpoint completeness; agreement of downstream manuscript conclusions requires the separate downstream checks. See the [figure and validation map](../MANUSCRIPT_REPRODUCTION.md).

## Timing and historical provenance

The [timing-pilot record](PILOT_VALIDATION.md) preserves the 20-path and 100-path selections, comparison commands and original SGE instructions. It is a timing/provenance reference, not the entry point for a new full run.

`simulate_single_gene.py`, `combine_rums_profile.py` and the SGE shell script are historical sources. Use the current portable runner and audited exporter for Sinabro 3.0.0. Earlier package imports and fixed cluster paths in historical scripts are not installation instructions for the current package.
