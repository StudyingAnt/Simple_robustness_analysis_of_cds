# Portable simulation and parallel validation

The Python entry points have no SeaWulf paths, module commands, or required scheduler environment variables. Cluster-specific settings belong in the Slurm wrappers.

- run_rums.py: simulation, optional input/signature subsets, deterministic sharding, multiprocessing and resume.
- validate_parallel.py: a bounded comparison of one worker versus multiple workers, using run_rums.py for both.
- pilot_selection.py: the same deterministic 12-transcript selection used by the earlier serial timing pilot.

Use [CONDA_SETUP.md](../CONDA_SETUP.md) to install Sinabro 3.0.0 and the recorded library versions. Run examples from the analysis repository root. Replace /PATH/TO paths with paths on your computer.

## Small reproducibility check without Slurm

The input is the complete preprocessed GENCODE v40 CDS-plus FASTA. The controller selects 12 transcripts, three per fixed length bin, then uses only SBS1, SBS2 and SBS13. The largest transcript is included. Internal-stop transcripts are excluded from this sample only; the preprocessing dataset is unchanged.

~~~bash
python 02_run_sinabro/validate_parallel.py \
  --input /PATH/TO/PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --out-dir /PATH/TO/PROJECT/runs/parallel_check_local \
  --sinabro-source ../sinabro \
  --workers 8 --n-traj 100 --seed 42
~~~

The output directory must be new, preventing a resumed run from being mistaken for a fresh timing measurement. Each mode executes 12 x 3 x 100 = 3,600 trajectories: 7,200 total, rather than the full human dataset.

Outputs include selected_transcripts.fa/csv, validation_manifest.json, serial/, parallel_8/, and comparison.json. Each mode uses the same input, seed and computation code. The comparison checks each accepted-step count and terminal code, not every intermediate DNA sequence. The earlier serial pilot separately audited intermediate sequences.

A passing comparison requires identical results, no incomplete paths, and observed overlap of at least two simulation cases. Worker PIDs and case start/end times are recorded in each mode's shards/*.cases.*.jsonl. comparison.json reports timing, simultaneous cases and elapsed-time speedup. It does not require eightfold speedup: process startup and long-transcript imbalance limit scaling. One-worker mode runs first, so timing may reflect cache/order effects; this is a correctness/concurrency smoke test, not a controlled performance study.

## Run on your own FASTA or an existing selected sample

The runner accepts CDS-plus FASTA records with unique transcript identifiers (the field before the first |). It does not require the full GENCODE count unless --expected-records is provided. CDS coordinates remain [1, len(sequence)-1): one flanking nucleotide at each end.

~~~bash
python 02_run_sinabro/run_rums.py \
  --input /PATH/TO/selected_transcripts.fa \
  --out-dir /PATH/TO/PROJECT/runs/custom_small_run \
  --workers 4 --n-traj 100 --seed 42 \
  --signatures SBS1 SBS2 SBS13
~~~

This example uses the normally installed package. --sinabro-source optionally selects a checkout. All 79 COSMIC v3.3.1 GRCh37 signatures are used if --signatures is omitted. The criterion is first protein-sequence change under the standard translation table.

Workstation defaults are one worker, one shard and shard index zero. The default trajectory count is the manuscript setting of 1,000; set --n-traj explicitly for small tests.

## Parallelism and sharding

One worker calculates one transcript/signature combination at a time; its independent trajectories run sequentially. Workers receive the next available combination. Each worker's NumPy/BLAS thread count is one. The default multiprocessing method is spawn with explicit initialization; no inherited fork-only state is required.

For multiple computers or scheduler array jobs, split the same input using --n-shards N and give each job a different --shard-index between 0 and N-1. Each machine runs its local pool; Python does not launch remote machines or require MPI.

The 100-shard full SeaWulf wrapper additionally passes --expected-records 73214. Its default remains at most eight concurrent jobs, each with one 96-CPU node. See [HPC_SUBMISSION.md](HPC_SUBMISSION.md) for the separate small SeaWulf job.

## Reproducibility, resume and platform scope

Checkpoints and run_manifest.json fingerprint the input, signature table, Python runner, package source, scientific settings and Python/NumPy/Biopython versions. Worker count, current shard index and process method are execution settings and do not change trajectory seeds. Changing shard count or selected signatures requires a new output directory.

Rerun the same command/output directory to skip validated completed transcripts. Unfinished transcripts are recomputed. A profile with any incomplete path receives NaN statistics; incomplete accepted-step counts are -1. No-overwrite checkpoint and lock protections remain.

The renamed runner uses identity rums-v3 and a new source fingerprint. It rejects checkpoints from earlier versions; preserve old outputs and use a fresh directory.

Python has no scheduler dependency. It uses spawn and standard-library POSIX/Windows file-lock branches. This version has been exercised on Linux/WSL; native macOS and Windows have not been verified. Output storage must support file locks and hard links for safe concurrent checkpoint publication. WSL provides the Linux workflow on Windows.

This checks the revised implementation's execution and reproducibility. It does not establish reproduction of every historical manuscript figure or complete the bat/R analysis stages.
