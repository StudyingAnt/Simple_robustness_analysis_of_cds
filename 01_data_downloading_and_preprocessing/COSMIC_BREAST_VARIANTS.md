# Breast-cancer observations mapped to the simulated BRCA1 transcript

This stage uses the two downloaded **COSMIC v101 / GRCh38** archives, the existing pinned GENCODE v40 CDS-plus FASTA and a small reference coordinate map included in Git. No COSMIC login, network download, GTF transfer or new simulation is needed on the compute node.

## Submit on SeaWulf

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/01_data_downloading_and_preprocessing/prepare_cosmic_breast_variants_seawulf.slurm
```

The job requests **one CPU, 4 GB, 30 minutes** (allocation limits, not measured runtime). It activates the existing `PROJECT/sinabro-revision` environment.

Default inputs:

```text
PROJECT/data/raw/cosmic/v101/GRCh38/Cosmic_MutantCensus_Tsv_v101_GRCh38.tar
PROJECT/data/raw/cosmic/v101/GRCh38/Cosmic_Classification_Tsv_v101_GRCh38.tar
PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa
```

Here `PROJECT` is the project directory from which `sbatch` is submitted. The existing first extraction job/output directory is not required. Both raw TARs and the FASTA are checked against their pinned hashes. If the user's `/data/raw/...` refers to the literal filesystem root rather than the project directory, override `COSMIC_RAW_DIR`:

```bash
COSMIC_RAW_DIR=/data/raw/cosmic/v101/GRCh38 \
  sbatch Simple_robustness_analysis_of_cds/01_data_downloading_and_preprocessing/prepare_cosmic_breast_variants_seawulf.slurm
```

Other overrides: `PROJECT`, `ANALYSIS_DIR`, `SINABRO_ENV`, `SINABRO_PYTHON`, `SINABRO_FASTA`, `COSMIC_BREAST_OUTPUT`.

## Selection and coordinate rules

1. Read the official compressed TSVs directly from the verified TARs. Retain exact `GENE_SYMBOL=BRCA1` observations.
2. Join `COSMIC_PHENOTYPE_ID` to Classification. Select `PRIMARY_SITE` equal to `breast`, ignoring case and surrounding whitespace. Do not infer tissue from a sample name or missing metadata. Conflicting tissue assignments fail. Missing tissue/phenotype linkage for a BRCA1 observation is audited and prevents a passing result.
3. Require one genomic position (`GENOME_START=GENOME_STOP`) and distinct single A/C/G/T reference/alternate alleles. Keep all observed substitution classes at this stage. Indels, multi-base/unknown alleles and positions outside the target CDS are counted as exclusions. Do not filter by missing somatic-status, pathogenicity or driver labels.
4. Map GRCh38 genomic positions to the **ENST00000357654.9 / BRCA1-203** target CDS through the frozen GENCODE v40 map. Include the terminal stop codon, matching the simulated 5,592-nt sequence. Genomic alleles are treated as forward-strand reference/alternative; complement them for this minus-strand target. Do not use source `MUTATION_CDS` or source transcript CDS numbering.
5. Verify the mapped reference base against the exact simulated CDS-plus sequence. Where a simple genomic SNV `HGVSG` annotation is present, cross-check its position and alleles against the genomic columns. Reference mismatches, contradictory annotations or invalid genomic positions prevent a passing result, with details retained in audit output.
6. Deduplicate accepted variants by target CDS position/reference/alternative. Preserve source-observation counts, unique COSMIC sample-ID counts and genomic COSV IDs. These counts do not establish independent patients or biological replicates.

The source report has **3,037 BRCA1 observations**, all annotated on ENST00000471181.7 / BRCA1-210 (5,655 coding nt). The target CDS is different. A fixed 63-nt shift is not used: mapping follows each genomic coding interval and strand.

## Coordinate-reference provenance and local verification

The small [coordinate map](brca1_gencode_v40_coordinates.json) is derived from the [GENCODE v40 comprehensive chromosome GTF](https://www.gencodegenes.org/human/release_40.html), combining target `CDS` and `stop_codon` features. GTF coordinates are 1-based; see the [GENCODE format specification](https://www.gencodegenes.org/pages/data_format.html). There are **22 contiguous coding segments**, totaling **5,592 nt**, ordered in transcript direction.

The official GTF gzip MD5 (`14a867b82917c8c3006838c3a5053a3e`) was matched against the [release MD5 manifest](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/MD5SUMS). Its SHA256 is `6254e1c52470d74c37f5a7969ee5b1b7debdd9d6bbd9bf722e65fc3d873f3104`.

An independent genomic FASTA was retrieved from the fixed GRCh38 chromosome accession **NC_000017.11**, bases **43,045,678–43,124,096**. Assembling the frozen coding intervals in transcript orientation reproduced the entire pinned target CDS exactly: SHA256 `2a8bfbf3a33c29bc232c0624f2e685fbe1ca2f114d7e452e1e09d87c597e1929`. The [independent validation record](brca1_gencode_v40_coordinate_validation.json) contains the public NCBI retrieval URL and genomic-sequence hash. Only coordinate/provenance JSON is included in Git; the downloaded reference files are not required on HPC.

To rebuild the map independently, download the GTF and the NCBI region FASTA using their recorded URLs, then run with new output paths:

```bash
python 01_data_downloading_and_preprocessing/build_brca1_coordinate_map.py \
  --gtf /PATH/TO/gencode.v40.annotation.gtf.gz \
  --genomic-region /PATH/TO/NC_000017.11_brca1_region.fa \
  --out-file /PATH/TO/NEW_MAP.json \
  --validation-out /PATH/TO/NEW_VALIDATION.json
```

Local synthetic tests cover both strands, coding-segment boundaries, terminal stop positions, alternative source-transcript coordinates, phenotype joins, duplicate sample observations, reference mismatches, conflicting annotations and existing-output preservation. Tests do not redistribute COSMIC observations. The real author TARs remain on HPC; the author subsequently supplied the successful full-data report below.

### Completed author-reported HPC preparation

The full-data report passed with no attention reasons at `2026-10-22T02:41:35.366043+00:00` (2026-10-21 in New York). Of 1,984,965 source rows, 3,037 were BRCA1. Classification selected **354 breast observations**, all source-annotated on ENST00000471181.7. Of these, **235 observations mapped as SNVs** to target ENST00000357654.9; they represent **174 distinct substitutions**. Exclusions: **74 non-SNV/unknown alleles** and **45 positions outside the target CDS**. The remaining 2,683 BRCA1 observations were from other primary sites. These are counts of catalog observations, not independent patients.

Both official archive hashes, the pinned FASTA, processing code (`2cb3fee027803cdf8e78a91ec42c6cc312febc52da67e7f42f1b9f8e908436c1`) and coordinate map (`39b7ecca50d24fafca787121f9ec9079da03c7e051a1390e826255077160e89e`) matched. Reported `brca1_breast_target_snvs.tsv` SHA256: `48bdae70c6a2417a4545fe4ae7341cf4abf106b81b0decdca2a6f1f0419cf521`; observation audit SHA256: `c53b73468ccca61c5f1a7d00230b18217897c763b56a4fd1571a50d88194faa5`.

No unresolved phenotype linkage, reference mismatch, genomic annotation conflict or invalid genomic coordinate was reported. The raw/processed observations remain on HPC and have not been inspected locally. The supplied report did not include its Slurm job ID/output-directory suffix. **Next: [match the SNVs to the existing graph](../07_manuscript_validation/BRCA1_OBSERVED_SNVS.md)**; the successful mapping job does not need resubmission.

## Outputs and result check

Default output:

```text
PROJECT/data/processed/cosmic_v101_grch38/brca1_breast_JOBID/
```

- **`cosmic_breast_variants_summary.json`**: pass/attention status, source and tissue counts, mapping/exclusion counts, distinct target variants and all input/reference/source hashes.
- `brca1_observation_audit.tsv`: all BRCA1 source rows with tissue assignment and mapping/exclusion status. Original columns and duplicate observations remain available.
- `brca1_breast_target_snvs.tsv`: one row per accepted target CDS substitution, with observation/sample-ID counts and COSV identifiers.

Share the summary JSON. `preparation_passed: true` requires accepted mapped SNVs and no unresolved BRCA1 phenotypes, reference mismatches, genomic annotation conflicts or invalid genomic coordinates. Biological exclusions such as non-breast observations, indels and positions outside the target CDS are reported separately. If attention is required, the summary/audit files are still saved and the Slurm job exits unsuccessfully to make that visible. Hash/schema errors stop before publishing processed outputs.

**This stage prepares the observed SNV set; it does not yet draw an observed-variant overlay.** A matched single substitution does not mean that an entire simulated multi-mutation DNA state was observed in a tumour. The historical observation-to-node rule remains unrecovered. Any regenerated annotation must state its matching definition. COSMIC observations alone do not identify driver status or establish that APOBEC caused a mutation. Existing simulated paths remain reusable.

## Portable execution

With the existing pinned environment active, from the analysis checkout:

```bash
python 01_data_downloading_and_preprocessing/prepare_cosmic_breast_variants.py \
  --raw-dir /PATH/TO/cosmic/v101/GRCh38 \
  --input /PATH/TO/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --out-dir /PATH/TO/NEW_PROCESSED_DIRECTORY
```

The Python implementation uses the standard library. Slurm configuration is isolated in the SeaWulf wrapper. Keep the restricted observations, audit rows and processed variant sets outside Git.
