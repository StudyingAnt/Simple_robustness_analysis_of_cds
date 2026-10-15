# SeaWulf submission and resume

The single file run_rums_seawulf.slurm contains the Slurm directives and Python runner. It requires the preprocessed GENCODE v40 FASTA, the Sinabro 3.0.0 source checkout, and a Python environment with NumPy and Biopython.

Default: short-96core, 100 array shards, at most 8 running shards, one 96-core node per shard, 4 hours per shard. Total concurrency is at most 8 nodes / 768 workers. All 73,214 transcripts and 79 COSMIC v3.3.1 GRCh37 profiles are processed with 1,000 independent trajectories per combination, seed 42, protein identity criterion, CDS start 1. Internal-stop records are retained as in the original preprocessing.

The limits follow the [SeaWulf queue guide](https://rci.stonybrook.edu/docs/queue-selection-guide). Actual allocation depends on scheduler availability and account access.

## Submit on the cluster

Log in through milan.seawulf.stonybrook.edu. First follow [Conda environment setup](../CONDA_SETUP.md), then transfer or regenerate the FASTA on the cluster. The package branch work/package-refactor-2026-10-14 contains Sinabro 3.0.0.

From the analysis directory:

```bash
conda activate sinabro-revision
mkdir -p logs
export SINABRO_PYTHON="$CONDA_PREFIX/bin/python"
export SINABRO_SOURCE=/absolute/path/to/sinabro
export SINABRO_FASTA=/absolute/path/to/gencode.v40.pc_transcripts.nopary.cdsplus.fa
export SINABRO_OUTPUT=/absolute/path/to/results/hpc_gencode_v40_n1000
sbatch 02_run_sinabro/run_rums_seawulf.slurm
```

The environment name matches the analysis repository's environment.yml; use your own name if you selected a different one. Recorded versions: Python 3.9.15, NumPy 1.21.6, Biopython 1.79. The script records actual versions and rejects mixing them within one output directory. It does not install packages. SINABRO_PYTHON selects the activated environment's absolute Python executable; the environment and all input/source paths must be accessible from compute nodes.

For a standalone copy, submit sbatch /path/to/run_rums_seawulf.slurm from a directory containing logs/. No other analysis Python script is required.

## Resume

After the previous submission ends, resubmit the same command with the same variables and output directory. Completed transcript checkpoints are validated and skipped. An unfinished transcript is recalculated; individual in-progress trajectory states are not checkpointed. Existing checkpoint files are not overwritten.

To retry only failed shards while retaining the original 100-shard assignment:

```bash
sbatch --array=3,17,52%8 02_run_sinabro/run_rums_seawulf.slurm
```

Per-shard locks reject simultaneous duplicates. A different input, profile table, seed, trajectory count, package/runner source, or library version requires a fresh SINABRO_OUTPUT directory.

## Outputs

- run_manifest.json: configuration, versions and input/code/profile hashes.
- transcripts/ENST*.npz: accepted-step matrix (1000, 79), termination codes, complete/incomplete counts, mean, SD, median, Monte Carlo SE and metadata.
- shards/NNN.done.json: computational completion of a shard.
- logs/: Slurm progress and errors.

Load checkpoints with numpy.load(path, allow_pickle=False). Signature order is in metadata. Terminal codes are 0: first nonsynonymous event, 1: no available mutation, 2: safety cap. Incomplete paths have accepted_steps=-1; a profile containing any incomplete path has NaN aggregate statistics. A finished shard means all cases were attempted, so inspect incomplete counts before downstream analysis.

Compressed arrays avoid billions of verbose CSV rows. Downstream manuscript scripts will require an explicit combination/export step for this format. This file submits the simulation stage; it does not generate manuscript figures.

Local verification covers Bash/Python syntax, checkpoint round trips, configuration mismatches, overwrite protection and incomplete-path handling. No HPC job or additional trajectory was executed during script creation.
