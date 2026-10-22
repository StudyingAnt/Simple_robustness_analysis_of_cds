# Running RUMS on your computer

Start with the runnable commands in [README.md](README.md). They use the installed Sinabro package, explicit input/output paths and ordinary Python processes.

## Input and parameters

`run_rums.py` accepts CDS-plus FASTA: one flanking base on each side, with coding sequence `sequence[1:-1]`. Identifiers must be unique before the first `|`. The human table exporter additionally requires the GENCODE transcript/gene IDs and names; the bundled synthetic example has this header format.

`--workers` is the number of local worker processes (default 1). `--n-traj` is paths per transcript/signature (default 1,000). `--seed 42` is the manuscript master seed. `--signatures SBS1 SBS2` selects a subset; omission uses all 79 fixed profiles. `--sinabro-source /PATH/TO/sinabro` optionally selects a package checkout instead of the installed package. Inspect all options with `python 02_run_sinabro/run_rums.py --help`.

One process handles one transcript/signature case at a time. Worker count does not alter the per-path random streams. The default process start method is `spawn`; each worker uses one numerical-library thread.

## Saved outputs and resume

Each simulation directory contains `run_manifest.json`, `transcripts/*.npz` and `shards/*.done.json`. Manifests fingerprint the input, signature table, scientific settings, package and runner source, and recorded runtime versions. Repeating exactly the same simulation command/output directory validates and skips completed transcripts; interrupted transcripts are recomputed. Changed scientific settings require a new directory. Audits, exports and figure generation require new output directories and do not resume simulations.

The runner also accepts `--n-shards N --shard-index K` to split a dataset into deterministic disjoint portions. The ordinary single-computer commands use the default one shard. The optional author execution archive records how the full calculations were distributed.

## Audit before export

`validate_results.py` checks the exact file set, configuration/input fingerprints, dimensions, terminal codes and recomputed per-path summary statistics. Inspect `validation_summary.json`, `transcript_audit.csv` and `file_set.json`. Success requires a complete passing audit; done markers alone are insufficient. [AGGREGATION.md](AGGREGATION.md) gives the full commands and table definitions.

The workflow has been exercised on Linux/WSL. Native Windows/macOS execution is not independently verified here. Checkpoint storage must support file locks and hard links. Numerical reproducibility is checked with the pinned environment and source fingerprints; identical PDF bytes across platforms are not promised.

## Author execution record

The authors used SeaWulf; account access is institution-specific. Their [Slurm settings and execution records](author_hpc/README.md), including pilot checks, are optional reference material. They are not prerequisites for the commands above.
