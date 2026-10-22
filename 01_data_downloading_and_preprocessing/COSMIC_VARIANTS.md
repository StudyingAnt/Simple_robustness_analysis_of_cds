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

The [official download page](https://cancer.sanger.ac.uk/cosmic/download/cosmic) provides **Scripted download** in the selected file's popup. The current site's published instructions use Basic authentication to request a temporary signed URL from:

```text
https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=...&bucket=...
```

Use this API rather than the older `/cosmic/file_download/...` examples. An authorized COSMIC account is required. Browser login is separate from command-line authentication. Never put passwords, Basic authorization strings or signed URLs in Git, shell arguments or Slurm scripts. The downloader prompts for email and a hidden password, keeps authentication in memory, and sends it only to the fixed COSMIC API. Storage requests carry no account Authorization header. Redirects are rejected; signed URLs are not printed or included in receipts.

On SeaWulf, run in an **interactive terminal**:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
conda activate /gpfs/projects/DamaghiGroup/projects/sinabro_pub/sinabro-revision

PROJECT=/gpfs/projects/DamaghiGroup/projects/sinabro_pub
ANALYSIS_DIR="${PROJECT}/Simple_robustness_analysis_of_cds"
python "${ANALYSIS_DIR}/01_data_downloading_and_preprocessing/download_cosmic_variants.py" \
  --out-dir "${PROJECT}/data/raw/cosmic/v101/GRCh38"
```

Enter account details only at the terminal prompts. No new installation is needed: the script uses the Python standard library (Python 3.9 or later). On other computers, use the local checkout and replace `--out-dir` with the desired absolute raw-data directory. The Python script has no SeaWulf paths or scheduler settings.

For an offline plan, without login or network:

```bash
python 01_data_downloading_and_preprocessing/download_cosmic_variants.py --plan
```

If authentication or archive access is unavailable, the program stops; it does not substitute another release. During preparation, the API route and public metadata were inspected, an unauthenticated request returned HTTP 400, and offline tests passed. **An authenticated end-to-end download has not yet been verified.**

## Browser and temporary-link alternatives

In the browser select **v101**, the product, **GRCh38**, then **Download in browser**. Save both original `.tar` archives in the raw directory above, then verify without account input:

```bash
python "${ANALYSIS_DIR}/01_data_downloading_and_preprocessing/download_cosmic_variants.py" \
  --out-dir "${PROJECT}/data/raw/cosmic/v101/GRCh38" --verify-only
```

Alternatively add `--signed-links` to the normal command to paste each archive's temporary URL at a hidden terminal prompt instead of entering account credentials. Obtain URLs through the official scripted-download instructions. The site states that links expire after one hour. Expected storage hosts, filenames and SHA256 values are checked. Do not publish temporary URLs.

## Outputs and next step

Archives stay compressed and unchanged outside the code checkout. Each successful run writes `download_verification_<id>.json` with the pinned release/build, official hashes, source-manifest hash and verification time. Existing matching archives are reused without login/network; existing mismatched files are preserved and cause an error. Interrupted or incorrect downloads are not published as final archives.

Next inspect TSV headers, select BRCA1 breast-cancer single-nucleotide substitutions, and check coordinates/reference alleles against **ENST00000357654.9**. Record provenance, exclusions and the annotation definition explicitly. Existing 1,000-path simulations do not need repeating. Keep downloaded COSMIC tables in raw data; this repository distributes retrieval code and public checksums.
