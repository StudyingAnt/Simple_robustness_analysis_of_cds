# Prepare the human reference

Activate the [pinned environment](../CONDA_SETUP.md), then run these commands from the analysis repository root. `PROJECT` is an absolute directory on your own computer, outside this source checkout. No HPC account is required.

## Download GENCODE v40

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project

bash 01_data_downloading_and_preprocessing/download_gencode.sh \
  --out-dir "$PROJECT/data/raw/gencode_v40" --release 40
```

Source: [GENCODE human release 40](https://www.gencodegenes.org/human/release_40.html), GRCh38.p13 protein-coding transcript FASTA. The script uses the official EBI HTTPS endpoint, validates gzip integrity and writes `gencode.v40.pc_transcripts.fa.download.json` with URL, hashes, sizes and record count. The expected source contains **108,573** transcripts. Use this fixed release for manuscript reproduction.

## Extract CDS-plus

```bash
python 01_data_downloading_and_preprocessing/extract_cds_plus.py \
  --input "$PROJECT/data/raw/gencode_v40/gencode.v40.pc_transcripts.fa" \
  --out-dir "$PROJECT/data/processed/gencode_v40" --release 40
```

Expected output: **73,214** records in `gencode.v40.pc_transcripts.nopary.cdsplus.fa`, plus `preprocessing_manifest.json` and rejected-record files. The coding sequence must have length divisible by three, start with ATG, end with TAA/TAG/TGA and retain one nucleotide on each side. PAR_Y transcripts are excluded. Accepted CDS is `sequence[1:-1]`; internal stop codons are retained and flagged rather than silently filtered.

GENCODE header CDS coordinates are 1-based inclusive. For bounds `start` and `end`, coding sequence is `sequence[start-1:end]` and CDS-plus is `sequence[start-2:end+1]`. Failed-rule masks are 1: length, 2: start, 4: stop, 8: PAR_Y; missing flanks are recorded separately. Scripts refuse to overwrite completed output directories.

Next: [human RUMS calculation](../02_run_sinabro/README.md). The independent supplementary BRCA1 observed-variant workflow uses [COSMIC download instructions](COSMIC_VARIANTS.md) and [breast-SNV coordinate mapping](COSMIC_BREAST_VARIANTS.md).

## Author execution

The authors ran preparation on Stony Brook University's SeaWulf HPC system. It requires institutional access; the terminal commands above do not. Original site-specific settings and observed-variant execution records are in [author_hpc/](author_hpc/README.md). Adapt these examples to your own site if you use an HPC system.
