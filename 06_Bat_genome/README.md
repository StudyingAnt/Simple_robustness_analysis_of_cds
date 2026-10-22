# Bat references, COSMIC controls and experimental APOBEC profiles

This stage is separate from the completed human COSMIC run. The manuscript uses Molossus molossus assembly HLmolMol2 and Myotis myotis assembly HLmyoMyo6 from the [Bat1K Pilot Project](https://bds.mpi-cbg.de/hillerlab/Bat1KPilotProject/). Keep these assemblies and annotation versions for the initial reproduction rather than silently replacing them with newer assemblies.

For the checked download -> CDS extraction -> bounded serial/parallel sequence, follow [HPC_REPRODUCTION.md](HPC_REPRODUCTION.md), one stage at a time. After the pilot passes, [COSMIC_PRODUCTION.md](COSMIC_PRODUCTION.md) describes full eight-signature execution, archive auditing and bat-specific table export. Local real-reference extraction and pilot evidence are summarized in [archived bat reference checks](../historical/validation/BAT_REFERENCE_CHECKS.md). The reference full COSMIC audit passed (job 2327882); export (job 2328086) passed with both frozen transcript populations and zero incomplete paths. The separate [six-enzyme experimental workflow](experimental_apobec/README.md) passed its full three-dataset HPC audit (job 2336448); the table export (job 2338741) provides all 18 dataset/profile summaries without resimulation.

## Download raw references independently

The portable download_bat_inputs.py uses only the Python standard library. bat1k_pilot_sources.json records official URLs, expected index sizes and checksums. Both genome hashes are now pinned to a completed, gzip-validated official-source download on 2026-10-17; they are not independent publisher-provided checksums. Downloads check sizes, configured hashes and gzip CRC, then publish a new raw-data directory with download_manifest.json. No CDS extraction or simulation is performed by the downloader. A failed download publishes no completed directory; rerun into a new destination after fixing the failure.

From the HPC project root:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/download_bat_inputs_seawulf.slurm
```

Default destination: PROJECT/data/raw/bat1k_pilot_2020, which must not already exist. Set SINABRO_BAT_RAW to another new directory if needed. The wrapper activates PROJECT/sinabro-revision, requests one shared CPU, 4 GB and a four-hour limit, and downloads README, both allAnnotatedIsoforms BED12 files, and both compressed genome FASTAs. The two genome files total 1,394,004,198 bytes compressed; annotation/documentation add approximately 17 MB. The limit is not a measured runtime. This job can run independently of human PCA/statistics work.

On another computer, from the analysis repository root:

```bash
python 06_Bat_genome/download_bat_inputs.py \
  --out-dir /PATH/TO/PROJECT/data/raw/bat1k_pilot_2020
```

For the small files alone, use --annotations-only and a separate new output directory. Genomes remain compressed in this stage. CDS extraction and scientific input validation are separate steps.

## Historical cohort and checked CDS extraction

The historical make_sinabro_input_*_all.py scripts concatenate every BED exon without clipping to thickStart/thickEnd. BED block definitions describe exons, while thickStart/thickEnd delimit the coding segment in these gene annotations; see [UCSC BED documentation](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). A transcript can therefore have untranslated exon sequence outside its CDS.

Inspection of the official annotation snapshots downloaded on 2026-10-17 found:

| Annotation | BED transcripts | Outer transcript span differs from coding span |
| --- | ---: | ---: |
| HLmolMol2 | 55,697 | 23,975 |
| HLmyoMyo6 | 47,391 | 23,665 |

However, every ID in the historical result tables (31,185 Molossus and 23,213 Myotis transcripts) matches a BED record with equal transcript/CDS outer bounds. No outer-boundary UTR discrepancy was found in these retained cohorts. The versioned lists in cohorts/ preserve those populations explicitly. extract_bat_cds.py checks the source/cohort hashes and coordinate agreement, extracts coding exon blocks and genomic flanks on both strands, and requires every frozen transcript to pass sequence validation. It does not admit additional isoforms from the larger annotation. The original make_sinabro_input_* files remain unchanged.

Local extraction from both actual genomes passed with the original transcript counts. Internal stops were found in 383 Molossus and 79 Myotis transcripts; these remain in the full FASTAs and are flagged in QC, matching the historical filtering policy. The bounded pilot omits them and records that exclusion. The original historical FASTAs have not been recovered, so this is not a claim of byte identity with archived original inputs. Transcript IDs do not constitute gene or orthologue mappings.

NNC APOBEC preferences use two upstream bases relative to C, whereas the COSMIC SBS context is centered on the mutated base. The [experimental runner](experimental_apobec/README.md) checks strand expansion, mutation-target bounds and available flanks explicitly. It preserves the original one-base CDS flanks and requires the full NNC context. Experimental profiles are handled separately from the COSMIC runner.

## Analysis and provenance

The two bat COSMIC runs contain 435,184,000 paths in total. Applying all six experimental enzyme profiles to human and both bats contains 765,672,000 paths. Full HPC audit logs for both campaigns report no invalid or incomplete paths. These counts describe completed simulation scope, not biological replicates or evidence of agreement with every manuscript claim.

Use [the manuscript analysis workflow](../07_manuscript_validation/README.md) to compare the audited exports, generate PCA/statistical tables and plot the species panels. The native-process panel compares human COSMIC SBS2 with Molossus APOBEC3A and Myotis APOBEC3A3; applying a single fixed enzyme profile to all three genomes is a different comparison. Both are reported explicitly.

The owner-authorized source workbook is included under `experimental_apobec/source_data/`. Preparation stages a verified copy into project `data/raw/`, preserves its motif-average values and normalizes per-site NNC preferences. The exact historical probability CSV was not recovered, so the generated profiles are labelled as a documented reconstruction. See [profile preparation and execution](experimental_apobec/README.md) for the conversion policy, profile names, checksums and commands.

The historical transcript populations are frozen. Their IDs alone do not provide gene or one-to-one orthologue mappings. The original transcript-level comparisons must not be interpreted as matched orthologue tests or as independent biological replicates. Raw experimental reads and the upstream assay-to-workbook pipeline have not been reprocessed by this workflow.
