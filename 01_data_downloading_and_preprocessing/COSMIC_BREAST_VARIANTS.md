# Prepare COSMIC breast SNVs for the BRCA1 graph

Obtain both pinned [COSMIC v101 GRCh38 archives](COSMIC_VARIANTS.md) and prepare the [GENCODE v40 CDS-plus FASTA](README.md). From the analysis repository root with the pinned environment active:

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
python 01_data_downloading_and_preprocessing/prepare_cosmic_breast_variants.py \
  --raw-dir "$PROJECT/data/raw/cosmic/v101/GRCh38" \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --out-dir "$PROJECT/data/processed/cosmic_v101_grch38/brca1_breast"
```

Choose a new output directory. This reads the verified TARs directly, joins tissue classifications and maps accepted genomic SNVs to **ENST00000357654.9 / BRCA1-203**. No scheduler, institutional account or prior BRCA1 source-extraction output is required. COSMIC archive retrieval itself requires the appropriate COSMIC account.

Inspect `cosmic_breast_variants_summary.json`: `preparation_passed=true`, no attention reasons, correct target transcript, and passing strand/reference checks. With the pinned inputs, reference counts are 354 breast BRCA1 observations, 235 accepted mapped observations and 174 distinct SNVs; 74 non-SNV/unknown-allele and 45 outside-CDS observations are excluded. Source observations are not independent patients.

Outputs include `brca1_observation_audit.tsv` (per-observation inclusion/exclusion audit), `brca1_breast_target_snvs.tsv` (deduplicated target substitutions) and the summary with source/code hashes. Keep raw and processed COSMIC observations outside Git. Next: [annotate the saved BRCA1 graph](../07_robustness_analysis/BRCA1_OBSERVED_SNVS.md).

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

Local synthetic tests cover both strands, coding-segment boundaries, terminal stop positions, alternative source-transcript coordinates, phenotype joins, duplicate sample observations, reference mismatches, conflicting annotations and existing-output preservation. Tests do not redistribute COSMIC observations. The source archives are retrieved separately by each authorized user and are not redistributed in Git.


## Author execution

The authors processed these inputs on SeaWulf. Their [actual submission settings and verified full-data record](author_hpc/COSMIC_BREAST_VARIANTS.md) are retained separately. SeaWulf is restricted to authorized institutional users and is not required by the command above. Adapt the optional HPC example to your own site's configuration.
