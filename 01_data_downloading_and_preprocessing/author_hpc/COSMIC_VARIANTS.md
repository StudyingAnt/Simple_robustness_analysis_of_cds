# COSMIC observations for the BRCA1 graph

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

Download the raw observations independently of the completed S12 simulation. Required products: **Census Genes Mutations** (coding mutations, including BRCA1) and **Classification** (tissue/histology metadata). Join them using `COSMIC_PHENOTYPE_ID` when selecting breast samples. The gene-only Cancer Gene Census table cannot replace mutation observations.

## Release and reference

Use **COSMIC v101, GRCh38** as a date-based candidate. The recovered `robustness_sbs2_blosum_20.h5` has filesystem modification date 2025-05-12; `lrho_sbs2.h5` has date 2025-05-19. These are modification dates, not verified creation dates. The BRCA1 graph artwork and notebook also have modification dates on 2025-05-19. [Official release notes](https://cancer.sanger.ac.uk/cosmic/release_notes?genome=38) place v101 on 2024-11-19 and v102 on 2025-05-21. Thus v101 is the nearest earlier release under this selection rule. **This does not prove which release the original analysis used:** older downloaded data may have remained in use, and copying can change timestamps.

GRCh38 matches the assembly underlying GENCODE v40. The simulation's `v3.3.1 / GRCh37` SBS label identifies a mutation probability profile, not the build for coordinates of observed variants.

Public metadata was checked on 2026-10-21 and frozen in [cosmic_v101_grch38.json](../cosmic_v101_grch38.json):

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

On SeaWulf, use an interactive **Slurm compute-node allocation** for downloads:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
srun --job-name=cosmic-download --partition=short-96core-shared \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2G --time=00:30:00 \
  --pty bash -lc '
    set -eo pipefail
    module load anaconda/3
    source "$(conda info --base)/etc/profile.d/conda.sh"
    PROJECT=/gpfs/projects/DamaghiGroup/projects/sinabro_pub
    conda activate "${PROJECT}/sinabro-revision"
    python "${PROJECT}/Simple_robustness_analysis_of_cds/01_data_downloading_and_preprocessing/download_cosmic_variants.py" \
      --out-dir "${PROJECT}/data/raw/cosmic/v101/GRCh38"
  '
```

Enter account details only at the terminal prompts. No new installation is needed: the script uses the Python standard library (Python 3.9 or later). On other computers, use the local checkout and replace `--out-dir` with the desired absolute raw-data directory. The Python script has no SeaWulf paths or scheduler settings.

For an offline plan, without login or network:

```bash
python 01_data_downloading_and_preprocessing/download_cosmic_variants.py --plan
```

If authentication or archive access is unavailable, the program stops rather than substituting another release. Automated tests check request construction and file handling without live account credentials; account authentication is exercised by the person retrieving the data. The reference archives were acquired and passed their official SHA256 checks.

## Browser and temporary-link alternatives

In the browser select **v101**, the product, **GRCh38**, then **Download in browser**. Save both original `.tar` archives in the raw directory above, then verify without account input:

```bash
PROJECT=/gpfs/projects/DamaghiGroup/projects/sinabro_pub
ANALYSIS_DIR="${PROJECT}/Simple_robustness_analysis_of_cds"
python "${ANALYSIS_DIR}/01_data_downloading_and_preprocessing/download_cosmic_variants.py" \
  --out-dir "${PROJECT}/data/raw/cosmic/v101/GRCh38" --verify-only
```

Alternatively add `--signed-links` to the normal command to paste each archive's temporary URL at a hidden terminal prompt instead of entering account credentials. Obtain URLs through the official scripted-download instructions. The site states that links expire after one hour. Expected storage hosts, filenames and SHA256 values are checked. Do not publish temporary URLs.

## Outputs and next step

Archives stay compressed and unchanged outside the code checkout. Each successful run writes `download_verification_<id>.json` with the pinned release/build, official hashes, source-manifest hash and verification time. Existing matching archives are reused without login/network; existing mismatched files are preserved and cause an error. Interrupted or incorrect downloads are not published as final archives.

### First step after receiving the MutantCensus TAR

Only the MutantCensus TAR is needed for the first source-extraction job. It checks the official SHA256, inspects TAR membership and TSV headers, and streams **every exact-BRCA1 row** into a smaller TSV. It retains duplicate sample observations, indels and other transcripts. Nothing is assigned to a simulated node at this stage. The raw TAR is not changed and no network or authentication is needed.

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/01_data_downloading_and_preprocessing/author_hpc/prepare_cosmic_brca1_seawulf.slurm
```

Default input:

```text
PROJECT/data/raw/cosmic/v101/GRCh38/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar
```

If the file was saved elsewhere, pass its existing absolute path; moving or unpacking it first is unnecessary:

```bash
COSMIC_MUTANTCENSUS_TAR=/ABSOLUTE/PATH/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar \
  sbatch Simple_robustness_analysis_of_cds/01_data_downloading_and_preprocessing/author_hpc/prepare_cosmic_brca1_seawulf.slurm
```

The job requests one CPU, 4 GB and 30 minutes (allocation limits). Outputs are placed in **`PROJECT/data/processed/cosmic_v101_grch38/brca1_source_JOBID/`**:

- `cosmic_brca1_preparation_summary.json`: hash verification, archive members, headers, source/BRCA1 row counts and transcript/mutation-type counts.
- `brca1_source_rows.tsv`: complete source columns for the retained BRCA1 observations.
- `transcript_counts.csv`: observation counts per source transcript.

Check `logs/cosmic-brca1-input-JOBID.out` and share the summary JSON. `preparation_passed: true` verifies this extraction; **`breast_filter_applied: false` is intentional** because Classification has not been used. Local tests use synthetic archives rather than redistributed COSMIC observations.

### Completed author-reported HPC extraction

On 2026-10-21 (New York), the author supplied a successful HPC report with UTC timestamp `2026-10-22T02:25:46.298301+00:00`. The official archive SHA256 matched. The reported processing-script SHA256 (`3adb18b3f623a53699bc32b5951ae256e90e7eb4dd45b56b5e065ff584d1ddcf`) and source-manifest SHA256 (`ba7900d9e37b6d03b8a411be15344f468037353aab70b713a108ac31e145de1e`) also match the local files.

- **1,984,965 source rows** scanned; **3,037 BRCA1 observations** retained.
- All BRCA1 observations use source transcript **ENST00000471181.7**.
- Reported `brca1_source_rows.tsv` SHA256: `4d72c7d9c0c6e2fce51885308adeae286a687e51dbaa5d3f5ffadc30d8e281c9`.
- Reported `transcript_counts.csv` SHA256: `4ea951b03f1a2b4ce6d396583b4d90a017f4b130be546e465ce0f5d552ee1a83`.
- No tissue, SNV or target-coordinate filter has been applied. Counts represent source observations, not distinct variants or independent patients.

The source observations remain on HPC; their content has not been independently inspected locally. The report omitted its Slurm job ID/output-directory suffix. Successful source extraction does not require resubmission.

**Coordinate check:** the local pinned GENCODE v40 CDS-plus FASTA contains source transcript ENST00000471181.7 (**BRCA1-210, 5,655 coding nt**) and simulation target ENST00000357654.9 (**BRCA1-203, 5,592 coding nt**). Their CDS hashes differ. Source `MUTATION_CDS` positions therefore must not be copied directly into target graph coordinates. Use the GRCh38 genomic position/alleles with the target's GENCODE v40 CDS annotation, account for transcript strand and verify each mapped reference allele against the target sequence. The target sequence and completed graph do not need to change for annotation.

Classification is required to identify breast samples through `COSMIC_PHENOTYPE_ID`. The author subsequently supplied both archives and a [passing full-data breast/coordinate-mapping report](../COSMIC_BREAST_VARIANTS.md): 354 breast observations, 235 accepted mapped observations and 174 distinct SNVs, with no attention reasons. The [observed-allele graph annotation](../../07_robustness_analysis/BRCA1_OBSERVED_SNVS.md) also passed with 19 single-substitution matches and an explicit caption rule. Existing extraction/mapping jobs do not need resubmission. This first source-extraction stage itself does not perform those later analyses.

Portable equivalent with the existing environment active:

```bash
python 01_data_downloading_and_preprocessing/prepare_cosmic_brca1.py \
  --archive /PATH/TO/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar \
  --out-dir /PATH/TO/NEW_PROCESSED_DIRECTORY
```

Then use the Classification archive to select breast samples and examine single-nucleotide substitutions/reference alleles against **ENST00000357654.9**. Record exclusions and matching semantics explicitly. Existing 1,000-path simulations do not need repeating. Keep raw and processed COSMIC observations outside Git; this repository distributes retrieval/processing code and public checksums.
