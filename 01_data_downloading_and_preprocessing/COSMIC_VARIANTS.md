# COSMIC observations for the BRCA1 graph

Download the raw observations independently of the completed S12 simulation. Required products: **Census Genes Mutations** (coding mutations, including BRCA1) and **Classification** (tissue/histology metadata). Join them using `COSMIC_PHENOTYPE_ID` when selecting breast samples. The gene-only Cancer Gene Census table cannot replace mutation observations.

## Release and reference

Use **COSMIC v101, GRCh38** as a date-based candidate. The recovered `robustness_sbs2_blosum_20.h5` has filesystem modification date 2025-05-12; `lrho_sbs2.h5` has date 2025-05-19. These are modification dates, not verified creation dates. The BRCA1 graph artwork and notebook also have modification dates on 2025-05-19. [Official release notes](https://cancer.sanger.ac.uk/cosmic/release_notes?genome=38) place v101 on 2024-11-19 and v102 on 2025-05-21. Thus v101 is the nearest earlier release under this selection rule. **This does not prove which release the original analysis used:** older downloaded data may have remained in use, and copying can change timestamps.

GRCh38 matches the assembly underlying GENCODE v40. The simulation's `v3.3.1 / GRCh37` SBS label identifies a mutation probability profile, not the build for coordinates of observed variants.

Public metadata was checked on 2026-10-21 and frozen in [cosmic_v101_grch38.json](cosmic_v101_grch38.json):

| Product | Archive | Official SHA256 |
| --- | --- | --- |
| [Census Genes Mutations](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/mutantcensus) | `Cosmic_MutantCensus_Tsv_v101_GRCh38.tar` (about 91 MB) | `9464b5a72f90747031fb2d13e4b42fca121f2f2dbe31643d2cff756d93e2c499` |
| [Classification](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/classification) | `Cosmic_Classification_Tsv_v101_GRCh38.tar` (about 164 kB) | `54f9e8d31a5cd768cff9c7762dc7e286881b0c9ea0bb05105f20a9522b188135` |

The original source-table filename and observation-to-graph matching rule remain unrecovered. These archives support a documented reconstruction; downloading alone does not reproduce the observed-variant overlay.

## Interactive scripted download

### Direct file-selection links

1. [Sign in to COSMIC](https://cancer.sanger.ac.uk/cosmic/auth/login).
2. [v101 Census Genes Mutations](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/mutantcensus): select **Cosmic_MutantCensus_Tsv_v101_GRCh38.tar**.
3. [v101 Classification](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/classification): select **Cosmic_Classification_Tsv_v101_GRCh38.tar**.

Click the selected filename to open the download popup. The `/download/cancer-mutation-census` landing page describes the separate Cancer Mutation Census product; use the direct COSMIC product links above for these inputs.

### Official command-line workflow

The official v101 MutantCensus scripted-download instructions describe three steps:

1. Encode the registered email and password, joined by a colon, as Base64. This is a reusable account credential, not encryption. The Python downloader performs this in memory after hidden password input.
2. Request a temporary download URL with the Basic authentication header. The exact public request URL for this archive is [the v101 MutantCensus scripted-download API](https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v101/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar&bucket=downloads):

   ```bash
   # Template only; prefer the Python downloader's hidden credential prompts.
   curl -H 'Authorization: Basic <BASE64_AUTHENTICATION_STRING>' \
     'https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v101/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar&bucket=downloads'
   ```

3. Download from the returned URL **without** sending the account authentication header:

   ```bash
   curl '<TEMPORARY_DOWNLOAD_URL_FROM_STEP_2>' \
     --output Cosmic_MutantCensus_Tsv_v101_GRCh38.tar
   ```

The popup states that the temporary URL is valid for one hour. The placeholders above are deliberately not real credentials or live signed URLs. For actual terminal use, the Python downloader below avoids putting account credentials into command arguments or history.

### Classification: official command-line example

For **Cosmic_Classification_Tsv_v101_GRCh38.tar**, start at the [v101 Classification file-selection page](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/classification) and select that exact GRCh38 archive.

1. Generate the Base64 authentication string from the registered email and password as described above. The same account authentication applies to both products; the downloader handles this in memory.
2. Request the Classification archive's temporary URL using [the official Classification scripted-download API](https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v101/Cosmic_Classification_Tsv_v101_GRCh38.tar&bucket=downloads):

   ```bash
   # Template only; prefer the Python downloader's hidden credential prompts.
   curl -H 'Authorization: Basic <BASE64_AUTHENTICATION_STRING>' \
     'https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v101/Cosmic_Classification_Tsv_v101_GRCh38.tar&bucket=downloads'
   ```

3. Use the returned URL to save the Classification archive, without an account authentication header:

   ```bash
   curl '<TEMPORARY_CLASSIFICATION_DOWNLOAD_URL_FROM_STEP_2>' \
     --output Cosmic_Classification_Tsv_v101_GRCh38.tar
   ```

The URL expires after one hour. Save this archive alongside MutantCensus under `PROJECT/data/raw/cosmic/v101/GRCh38/`. Its official SHA256 is `54f9e8d31a5cd768cff9c7762dc7e286881b0c9ea0bb05105f20a9522b188135`. The popup's Actionability example illustrates URL syntax for another product; use the Classification path and output filename shown here for this analysis.

The existing Python downloader already includes both products. With a valid MutantCensus archive in the output directory, it verifies and reuses that file, then downloads only the missing Classification archive. After both files are present, `--verify-only` checks them without account authentication.

### API configuration and portable downloader

The official request uses `bucket=downloads`, the API alias. It is distinct from the physical storage bucket `cosmic-downloads-production`; both are recorded separately in the pinned manifest.

The [official download page](https://cancer.sanger.ac.uk/cosmic/download/cosmic) provides **Scripted download** in the selected file's popup. The current site's published instructions use Basic authentication to request a temporary signed URL from:

```text
https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=...&bucket=...
```

Use this API rather than the older `/cosmic/file_download/...` examples. An authorized COSMIC account is required. Browser login is separate from command-line authentication. Never put passwords, Basic authorization strings or signed URLs in Git, shell arguments or Slurm scripts. The downloader prompts for email and a hidden password, keeps authentication in memory, and sends it only to the fixed COSMIC API. Storage requests carry no account Authorization header. Redirects are rejected; signed URLs are not printed or included in receipts.

Run from the analysis repository root with your environment active:

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
python 01_data_downloading_and_preprocessing/download_cosmic_variants.py \
  --out-dir "$PROJECT/data/raw/cosmic/v101/GRCh38"
```

Enter account details only at the terminal prompts. The standard-library downloader needs Python 3.9 or later. To inspect pinned metadata without network or login, add `--plan` instead of `--out-dir`. To paste temporary links at hidden prompts, add `--signed-links` to the normal command.

## If you downloaded the TAR files in your browser

Save both unchanged archives under your chosen raw directory and verify them:

```bash
python 01_data_downloading_and_preprocessing/download_cosmic_variants.py \
  --out-dir "$PROJECT/data/raw/cosmic/v101/GRCh38" --verify-only
```

Each successful download/verification writes a receipt with pinned build/release, hashes and verification time. Existing matching files are reused; mismatching files are preserved and cause failure. Raw COSMIC observations and temporary credentials/URLs stay outside Git.

## Optional BRCA1-only source extraction

This streams every exact-BRCA1 observation from MutantCensus without a tissue or SNV filter. It is an inspection step; the full breast-coordinate preparation below can read both TARs directly.

```bash
python 01_data_downloading_and_preprocessing/prepare_cosmic_brca1.py \
  --archive "$PROJECT/data/raw/cosmic/v101/GRCh38/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar" \
  --out-dir "$PROJECT/data/processed/cosmic_v101_grch38/brca1_source"
```

Next: [select breast observations and map to the graph transcript](COSMIC_BREAST_VARIANTS.md). Source transcript CDS positions must not be copied into the target transcript without genomic mapping.

## Author execution

The authors downloaded and processed these inputs for calculations on SeaWulf. SeaWulf requires institutional access; COSMIC access separately requires an authorized COSMIC account. The [actual site-specific commands and verified source-extraction record](author_hpc/COSMIC_VARIANTS.md) are retained for provenance. Adapt those examples to your own HPC site if applicable; the portable commands above need no SeaWulf account.
