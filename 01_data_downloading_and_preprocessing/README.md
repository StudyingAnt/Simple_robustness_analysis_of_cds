# GENCODE v40 download and CDS-plus preprocessing

Run these commands from the analysis repository root. Downloaded and generated data are excluded from Git. Each stage writes a provenance manifest and refuses to overwrite existing outputs.

## 1. Download and unpack the transcript FASTA

```bash
bash 01_data_downloading_and_preprocessing/download_gencode.sh \
  --out-dir data/raw/gencode_v40 --release 40
```

Source: [GENCODE human release 40](https://www.gencodegenes.org/human/release_40.html), GRCh38.p13 protein-coding transcript FASTA. The script uses the official EBI HTTPS endpoint, checks gzip integrity, unpacks the file, and records hashes, sizes, record count, URL, UTC timestamps, and source-script hash in `gencode.v40.pc_transcripts.fa.download.json`.

The expected historical input contains 108,573 transcript records. This is a fixed historical release; do not substitute the latest GENCODE data when reproducing the manuscript.

## 2. Extract CDS-plus

Use Biopython 1.79 for comparison with the recorded environment. The full historical environment listing records Python 3.9.15, whereas this directory previously documented 3.9.16; retain that discrepancy in provenance. The current step is run with Python 3.9.15 and Biopython 1.79.

```bash
python 01_data_downloading_and_preprocessing/extract_cds_plus.py \
  --input data/raw/gencode_v40/gencode.v40.pc_transcripts.fa \
  --out-dir data/processed/gencode_v40 --release 40
```

The original selection rules are retained:

- CDS length is divisible by three.
- CDS begins with ATG and ends with TAA, TAG, or TGA.
- Records marked PAR_Y are excluded.
- One transcript nucleotide is present on each side of the CDS.

GENCODE header coordinates are 1-based inclusive. For CDS bounds `start` and `end`, Python slices are `sequence[start-1:end]` for CDS and `sequence[start-2:end+1]` for CDS-plus. The header is parsed by its CDS field rather than by a fixed field position.

The expected historical output is **73,214** CDS-plus records in `gencode.v40.pc_transcripts.nopary.cdsplus.fa`. Original filter failures are written to the fifteen `nopass.parity*.fa` files. Their bit masks are 1: length, 2: start, 4: stop, 8: PAR_Y; combined masks represent multiple failed rules. Otherwise-valid records lacking a flank are written separately to `nopass.missing_flank.fa`; the old script omitted these records without a separate file.

`preprocessing_manifest.json` records all output hashes, record counts, the input hash, environment versions, and rejection accounting. No new filter for ambiguous bases or internal stop codons is silently introduced.

## 3. Inspect before simulation

Verify the input/output counts, output hash against the historical reference when available, and the rejection accounting. Check that the accepted coding slice is `sequence[1:-1]`. The next stage is a small single-sequence simulation, handled separately from preprocessing.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests cover coding-coordinate boundaries, historical filter masks, both missing-flank cases, record accounting, invalid annotations, and preservation of existing output. The preserved original scripts remain on the backup branch and are used as an independent comparison for the fresh full dataset.
