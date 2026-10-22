# Bat CDS preparation and bounded COSMIC validation

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

This workflow starts from the raw references on an HPC project with no bat files.
It preserves the historical transcript cohorts and leaves the original scripts
unchanged. It does not submit full bat simulations or the experimental APOBEC,
CGC/BLOSUM or BRCA1 analyses. All Python scripts also run on other computers;
only the Slurm wrappers are specific to SeaWulf.

## 0. Update code and check prerequisites

Use the pinned environment and Sinabro 3.0.0 package checkout. The
environment setup is described in ../CONDA_SETUP.md. Scientific versions remain
Python 3.9.15, NumPy 1.21.6 and Biopython 1.79. No extra dependency is required.

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds switch work/reproduce-manuscript-2026-10-14
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
```

Submit every job below from this project root. The wrappers load anaconda/3,
activate PROJECT/sinabro-revision and obtain PROJECT from SLURM_SUBMIT_DIR unless
explicitly overridden. They do not install or modify the package/environment.
Do not chain job submissions until you have checked the preceding job's report.

## 1. Download raw references

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/download_bat_inputs_seawulf.slurm
```

Destination: data/raw/bat1k_pilot_2020/. The job retrieves the two compressed
assemblies, both BED annotations and the source README (about 1.4 GB compressed).
Success: job exits successfully and download_manifest.json has
download_passed=true, annotations_only=false and five outputs. Sizes, known
annotation hashes, pinned genome hashes and gzip integrity are checked/recorded.
Genome hashes were measured from the completed local official-source download on
2026-10-17 and pinned for subsequent runs; they are not publisher-provided checksums.

Existing destinations are never overwritten. If this stage has already passed,
use those files and skip downloading. An annotations-only directory is insufficient.
For another destination, set SINABRO_BAT_RAW consistently in steps 1 and 2.

## 2. Extract and validate both historical CDS cohorts

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/extract_bat_cds_seawulf.slurm
```

Destination: data/processed/bat1k_historical_v1/{molMol2,myoMyo6}/. The extractor
reads compressed genomes one contig at a time; no whole uncompressed genome copy
is required. It validates the raw-file sizes/hashes and any download manifest,
checks the versioned cohort lists, then uses BED12 exon/CDS intersections and
transcript-oriented flanks. Negative-strand transcripts are reverse complemented
including both flanks. Output bases are uppercased.

Each species directory contains:

- extraction_manifest.json: counts, input/output hashes, policies and pass status.
- transcript_qc.csv: coordinates, strand, exon count, CDS length, translated length,
  internal-stop count and errors for every historical transcript.
- SPECIES.cdsplus.fa: sorted transcript IDs and CDS plus one genomic base per end;
  published only when every historical transcript passes required checks.

Success requires extraction_passed=true, failed_transcripts=0 and exactly 31,185
Molossus / 23,213 Myotis records. Checks include cohort identity, CDS boundaries,
complete flanks, A/C/G/T bases, frame, ATG start and terminal stop. Internal stops
are retained and flagged, matching the historical filtering policy; their presence
does not by itself fail extraction. These are transcript counts, not gene counts.
The full cohort is not silently reduced if any required check fails.

The original scripts did not explicitly clip UTRs. Every historical retained ID
has equal transcript/CDS outer bounds in the frozen BED, so that difference does
not change these retained intervals. The new workflow does not add new isoforms.
Matching counts and passing new checks do not prove byte identity with unavailable
historical FASTAs. Failed sequence checks publish QC and a false manifest, without
a final FASTA. After investigating a failure, choose a new processed destination
with SINABRO_BAT_PROCESSED; successful earlier directories are preserved.

An independent local reconstruction of the old filters accepted two extra Myotis
IDs, mRNA9113 and mRNA9159, both with N-containing CDSs. They are absent from the
historical result tables and are not added here. Their historical exclusion reason
is unknown. All frozen-cohort sequences match the legacy extraction after uppercase
normalization; the Molossus old script left some flank letters lowercase.

## 3. Run a small serial/parallel check

Only after both extraction manifests pass:

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/validate_bat_pilot_seawulf.slurm
```

Destination: runs/bat_cosmic_pilot_JOBID/. For each species, select 12 evenly
spaced CDS-length quantile representatives, run SBS2, SBS4 and SBS13 with 100
trajectories per combination, and compare one worker with eight workers. Total:
7,200 trajectories per mode across both species, 14,400 including both modes.
An absolute max_steps=10,000 safety cap is used; a capped path is incomplete and
fails this validation. Internal-stop records are excluded only from this bounded
pilot, with exclusion counts recorded; they remain in the full input FASTA.

Success: pilot_report.json has validation_passed=true. For both species,
results_equal, all_paths_complete and observed_parallel_execution must be true;
both serial/parallel archive audits must pass. Counts and terminal reasons are
compared for every trajectory; archive integrity and saved summary statistics
are checked separately. Timings, observed worker PIDs/concurrency, selected IDs
and lengths, hashes and configuration are saved. Speedup is reported, but no
particular speedup factor is required. This demonstrates execution and consistency,
not scientific equivalence to the historical random paths or a whole-genome timing
estimate. No new experimental-enzyme profile is used in this test.

## Portable local commands

From the analysis repository root, with the pinned environment active:

```bash
python 06_Bat_genome/download_bat_inputs.py --out-dir /PATH/TO/PROJECT/data/raw/bat1k_pilot_2020
python 06_Bat_genome/extract_bat_cds.py --raw-dir /PATH/TO/PROJECT/data/raw/bat1k_pilot_2020 --out-dir /PATH/TO/PROJECT/data/processed/bat1k_historical_v1/molMol2 --species molMol2
python 06_Bat_genome/extract_bat_cds.py --raw-dir /PATH/TO/PROJECT/data/raw/bat1k_pilot_2020 --out-dir /PATH/TO/PROJECT/data/processed/bat1k_historical_v1/myoMyo6 --species myoMyo6
python 06_Bat_genome/author_hpc/validate_bat_pilot.py --processed-root /PATH/TO/PROJECT/data/processed/bat1k_historical_v1 --out-dir /PATH/TO/PROJECT/runs/bat_pilot_local --sinabro-source /PATH/TO/PROJECT/sinabro --workers 4 --records 12 --n-traj 20
```

Choose new output directories for a repeat. A failed or existing output is not
automatically deleted. After this pilot passes, continue with
[COSMIC_PRODUCTION.md](COSMIC_EXECUTION.md) for the full eight-signature array,
archive audit and bat-specific table export. The six experimental enzyme profiles
remain a separate workflow. aggregate_bat_results.py preserves bat IDs and reuses
the audited checkpoint reader; the GENCODE-only exporter must not be used directly.
