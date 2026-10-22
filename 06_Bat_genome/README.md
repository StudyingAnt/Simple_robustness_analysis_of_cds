# Bat CDS preparation and COSMIC calculations

Run from the analysis repository root with the [pinned environment](../CONDA_SETUP.md) active. The commands below run on your own computer without a scheduler or institutional account. Use an absolute `PROJECT` directory outside the source checkout and set `WORKERS` to your available CPU processes (1 is valid).

## 1. Download the original assemblies and annotations

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
WORKERS=1

python 06_Bat_genome/download_bat_inputs.py \
  --out-dir "$PROJECT/data/raw/bat1k_pilot_2020"
```

The reference is the [Bat1K Pilot Project](https://bds.mpi-cbg.de/hillerlab/Bat1KPilotProject/): HLmolMol2 (*Molossus molossus*) and HLmyoMyo6 (*Myotis myotis*). Here “Pilot” is part of the source project's name, not a required preliminary simulation. Two compressed genomes total about 1.4 GB; annotations add about 17 MB. The downloader verifies the pinned source hashes, sizes and gzip integrity and writes `download_manifest.json`.

## 2. Extract the frozen CDS-plus cohorts

```bash
for SPECIES in molMol2 myoMyo6; do
  python 06_Bat_genome/extract_bat_cds.py \
    --raw-dir "$PROJECT/data/raw/bat1k_pilot_2020" \
    --out-dir "$PROJECT/data/processed/bat1k_historical_v1/$SPECIES" \
    --species "$SPECIES"
done
```

Expected populations: **31,185** Molossus and **23,213** Myotis transcripts, as fixed by [cohorts/](cohorts/README.md). Outputs per species include `SPECIES.cdsplus.fa`, `transcript_qc.csv` and `extraction_manifest.json`. Every frozen transcript must pass extraction; the program does not select new isoforms. Internal-stop transcripts (383 Molossus, 79 Myotis) remain in the full inputs and are flagged.

## 3. Run the eight COSMIC signatures

```bash
for SPECIES in molMol2 myoMyo6; do
  python 02_run_sinabro/run_rums.py \
    --input "$PROJECT/data/processed/bat1k_historical_v1/$SPECIES/$SPECIES.cdsplus.fa" \
    --out-dir "$PROJECT/runs/bat_cosmic_n1000_seed42/$SPECIES" \
    --signatures SBS1 SBS2 SBS3 SBS4 SBS6 SBS9 SBS13 SBS18 \
    --n-traj 1000 --seed 42 --max-steps 10000 --workers "$WORKERS"
done
```

The complete calculation is **435,184,000 paths** and may take substantial time on a personal computer. This is the manuscript configuration, not a short example. The stopping criterion is the first protein-sequence change. The fixed COSMIC v3.3.1 GRCh37 label identifies the profile table; bat coordinates remain those of the original Bat1K assemblies. These direct commands do not require a pilot report.

## 4. Audit and export

```bash
python 06_Bat_genome/process_bat_results.py \
  --stage audit --processed-root "$PROJECT/data/processed/bat1k_historical_v1" \
  --run-root "$PROJECT/runs/bat_cosmic_n1000_seed42" \
  --out-root "$PROJECT/runs/bat_cosmic_audit" --n-traj 1000 --workers "$WORKERS"

python 06_Bat_genome/process_bat_results.py \
  --stage export --processed-root "$PROJECT/data/processed/bat1k_historical_v1" \
  --run-root "$PROJECT/runs/bat_cosmic_n1000_seed42" \
  --audit-root "$PROJECT/runs/bat_cosmic_audit" \
  --out-root "$PROJECT/runs/bat_cosmic_export" --n-traj 1000 --workers "$WORKERS"
```

Read `audit_manifest.json` and `export_manifest.json`: both should have `stage_passed=true`, with both full populations and no incomplete paths. Each species export contains signature summaries and transcript matrices. Then prepare the [six experimental enzyme profiles](experimental_apobec/README.md) and run the [cross-species manuscript analyses](../07_manuscript_validation/README.md).

## Method and interpretation

Extraction checks coding boundaries, exon order, strand, source hashes and frozen IDs. Historical exon-concatenation scripts did not clip every BED exon to its coding bounds, but all retained cohort IDs in these snapshots have equal transcript/CDS outer bounds. The current extractor checks this agreement and genomic flanks explicitly. Recovered cohort IDs do not establish byte identity to missing historical FASTAs.

All retained isoforms receive equal transcript weight. These are not one-to-one orthologue or gene-collapsed comparisons, and transcripts are not independent biological replicates. Experimental NNC profiles are a separate mutation model from the centered COSMIC SBS context.

## Author execution

The authors ran the full calculations on Stony Brook University's SeaWulf HPC system, which requires institutional access. Their actual Slurm files, completed-run records and pilot checks are in [author_hpc/](author_hpc/README.md). If you use your own HPC system, adapt these examples to your site's settings. SeaWulf access and its pilot reports are not prerequisites for the direct Python workflow above.
