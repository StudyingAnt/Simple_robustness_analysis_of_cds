# SeaWulf submission and resume

The simulation code is in [run_rums_seawulf.py](run_rums_seawulf.py). The [Slurm wrapper](run_rums_seawulf.slurm) requests resources, activates Conda, sets paths and thread limits, and calls the Python file with named arguments. The Python runner handles trajectory generation, node-local multiprocessing, deterministic shard assignment, and checkpoints.

Default: short-96core, 100 array shards, at most 8 running shards, one 96-core node per shard, 4 hours per shard. Total concurrency is at most 8 nodes / 768 workers. All 73,214 transcripts and 79 COSMIC v3.3.1 GRCh37 profiles are processed with 1,000 independent trajectories per combination, seed 42, protein identity criterion, CDS start 1. Internal-stop records are retained as in the original preprocessing.

The limits follow the [SeaWulf queue guide](https://rci.stonybrook.edu/docs/queue-selection-guide). Actual allocation depends on scheduler availability and account access.

**This is the full-workload runner.** This separation does not add a bounded pilot mode or establish HPC parallel performance. Selecting one array index still runs that shard of the full dataset across all 79 signatures.

## Project layout

Keep the code repositories separate from project data and execution outputs:

```text
PROJECT/
├── sinabro/
├── Simple_robustness_analysis_of_cds/
├── sinabro-revision/
├── data/
│   ├── raw/gencode_v40/
│   └── processed/gencode_v40/
├── runs/
├── results/
├── logs/
└── slurm/
```

The wrapper defaults to this layout. Inputs are read from data/processed; per-run checkpoints belong in runs. Reserve results for final combined tables and figures. The Conda environment path can be overridden if you store environments under envs/ or elsewhere.

## Submit on the cluster

Log in through milan.seawulf.stonybrook.edu. Follow [Conda setup](../CONDA_SETUP.md), and complete the download and preprocessing stages before submission. The package branch work/package-refactor-2026-10-14 contains Sinabro 3.0.0.

Submit from the **HPC project root**, not the analysis repository root. For the current project:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/02_run_sinabro/run_rums_seawulf.slurm
```

The wrapper loads anaconda/3 and activates the revision environment inside each array job. Adjust the module name to the one available on your cluster if necessary. It does not install packages. Recorded core versions: Python 3.9.15, NumPy 1.21.6, Biopython 1.79.

Optional exported variables:

| Variable | Default |
|---|---|
| PROJECT | Slurm submission directory |
| ANALYSIS_DIR | PROJECT/Simple_robustness_analysis_of_cds |
| SINABRO_ENV | PROJECT/sinabro-revision |
| SINABRO_SOURCE | PROJECT/sinabro |
| SINABRO_FASTA | PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa |
| SINABRO_OUTPUT | PROJECT/runs/hpc_gencode_v40_n1000_seed42 |
| SINABRO_PYTHON | Activated environment's bin/python |
| SINABRO_N_TRAJ | 1000 |
| SINABRO_SEED | 42 |
| SINABRO_MAX_STEPS | 10000 |

For a custom environment or run directory, export the corresponding variables before sbatch. Use absolute paths visible from compute nodes. The default output directory name incorporates the actual trajectory count and seed. Slurm log paths are relative to the submission directory; create logs/ there before submission.

You may copy the Slurm wrapper into PROJECT/slurm/ for site-specific changes; it still calls the Python file in ANALYSIS_DIR/02_run_sinabro/. Keep that checkout fixed while jobs run.

To inspect Python arguments without starting a simulation:

```bash
python Simple_robustness_analysis_of_cds/02_run_sinabro/run_rums_seawulf.py --help
```

## Resume

After the previous submission ends, resubmit the same command with the same variables and output directory. Completed transcript checkpoints are validated and skipped. An unfinished transcript is recalculated; individual in-progress trajectory states are not checkpointed. Existing checkpoint files are not overwritten.

To retry only failed shards while retaining the original 100-shard assignment:

```bash
sbatch --array=3,17,52%8 Simple_robustness_analysis_of_cds/02_run_sinabro/run_rums_seawulf.slurm
```

Per-shard locks reject simultaneous duplicates. A different input, profile table, seed, trajectory count, package/runner source, or library version requires a fresh SINABRO_OUTPUT directory. Worker count and the current shard index do not change the shared scientific configuration.

The separated runner records the SHA256 of the actual Python file and uses runner identity seawulf-rums-v2. Checkpoints from the previous embedded-Python runner have a different configuration fingerprint and are intentionally rejected. Preserve those outputs and choose a fresh directory; do not remove the configuration checks.

## Outputs

- run_manifest.json: configuration, versions and input/code/profile hashes.
- transcripts/ENST*.npz: accepted-step matrix (1000, 79), termination codes, complete/incomplete counts, mean, SD, median, Monte Carlo SE and metadata.
- shards/NNN.done.json: computational completion of a shard.
- logs/: Slurm progress and errors, outside the per-run output directory.

Load checkpoints with numpy.load(path, allow_pickle=False). Signature order is in metadata. Terminal codes are 0: first nonsynonymous event, 1: no available mutation, 2: safety cap. Incomplete paths have accepted_steps=-1; a profile containing any incomplete path has NaN aggregate statistics. A finished shard means all cases were attempted, so inspect incomplete counts before downstream analysis.

Compressed arrays avoid billions of verbose CSV rows. Downstream manuscript scripts require an explicit combination/export step for this format. The runner handles the simulation stage; it does not generate manuscript figures.

Local tests cover Bash/Python syntax, named arguments, checkpoint round trips, configuration mismatches, overwrite protection and incomplete-path handling. No HPC job or additional trajectory was executed during this separation.
