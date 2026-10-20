# Bat input preparation and remaining simulation scope

This stage is separate from the completed human COSMIC run. The manuscript uses Molossus molossus assembly HLmolMol2 and Myotis myotis assembly HLmyoMyo6 from the [Bat1K Pilot Project](https://bds.mpi-cbg.de/hillerlab/Bat1KPilotProject/). Keep these assemblies and annotation versions for the initial reproduction rather than silently replacing them with newer assemblies.

For the checked download -> CDS extraction -> bounded serial/parallel sequence, follow [HPC_REPRODUCTION.md](HPC_REPRODUCTION.md), one stage at a time. After the pilot passes, [COSMIC_PRODUCTION.md](COSMIC_PRODUCTION.md) describes full eight-signature execution, archive auditing and bat-specific table export. Local real-reference extraction and pilot evidence are summarized in [LOCAL_VALIDATION.md](LOCAL_VALIDATION.md). The user reported successful full COSMIC audit (job 2327882) and export (job 2328086), with both frozen transcript populations and zero incomplete paths. The separate [six-enzyme experimental workflow](experimental_apobec/README.md) passed its full three-dataset HPC audit (job 2336448); follow its export stage to obtain numerical tables without resimulation.

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

## Remaining long trajectory calculations

Under a descriptive revision that retains the bat results but withdraws unsupported adaptive claims:

| Work item | Additional trajectory simulation? | Scope/dependency |
| --- | --- | --- |
| Two bat genomes under the manuscript COSMIC controls | Completed; no repeat needed | Eight signatures: SBS1, SBS2, SBS3, SBS4, SBS6, SBS9, SBS13, SBS18. The human values for these are already available from the 79-signature run. User-reported full bat audit and table export passed. |
| Human plus both bats under six bat enzyme profiles | Completed; no repeat needed | The same six profiles on all three species, including own/foreign enzymes and human controls. User-reported full audit 2336448 passed 765,672,000 paths. Export and historical numerical comparison follow without resimulation. Exact historical probability CSVs remain unavailable; the workbook reconstruction policy is explicit. Human SBS2 is not a replacement for these profiles. |
| Monte Carlo sample-size sensitivity, SD/SE/zero fraction/CI | Usually no | Reuse saved per-path counts and N=100/250/500/1000 prefixes or resampling. A small independent-seed confirmation is possible; it does not require another whole-human run. |
| Weighted/unweighted one-step baselines, GC3/codon covariates, PCA, descriptive/inferential statistics, CGC comparison | No full RUMS rerun | Use sequence-derived features and existing matrices. Implement/check the revised statistical design separately. |
| BLOSUM criterion sensitivity | Limited subset if retained/needed | Reproduce the supporting sensitivity analysis on prespecified genes/signatures. Full-genome BLOSUM is not necessary solely to present a conservative sequence-identity criterion; long alignments can be expensive. |
| Synonymous-codon null | Conditional | Needed to support retained adaptation claims or if specifically required after reframing. Define a bounded prespecified design first rather than launching a whole-genome null by default. |
| Orthologue-matched species-by-process comparison | Conditional | Required to support stronger coevolution claims; define one-to-one orthologues/one transcript per gene first. Existing per-transcript results can be subsetted if the exact sequences are retained; new sequence representatives require only the missing cases. |

The frozen cohorts contain 31,185 Molossus and 23,213 Myotis transcripts, all successfully extracted locally. Two bats x eight signatures x 1,000 paths would be 435,184,000 paths. Six profiles on 73,214 human + 31,185 Molossus + 23,213 Myotis transcripts would be 765,672,000 paths. Their sum is approximately 21% of the completed human 79-signature path count. These are path-count illustrations, not runtime estimates; NNC preference, sequence length and path length affect runtime.

## Independently preparable tasks

1. Download the raw two-assembly references using the job above while human downstream analyses proceed.
2. Review the original six-enzyme workbook, map source IDs/species, preserve raw values and define normalization. The owner-authorized original workbook is included in experimental_apobec/source_data/; the preparation job stages a verified copy into project data/raw. Generated experimental profiles remain outside Git. The source workbook contains motif-average counts/percentages; they must not be silently treated as summed raw mutation counts. Source metadata, background/replicate information and the original conversion code determine which uncertainty claims can be made.
3. Follow HPC_REPRODUCTION.md to download, extract and check the small bat COSMIC pilot before requesting expensive jobs.
4. Prepare GC3/codon/one-step baselines and human PCA/statistical scripts against the completed human export, independently of the bat download/profile work.
5. Use the existing local CGC v99 GRCh38 snapshot for initial annotation/mapping checks rather than automatically changing releases or redownloading it. The historical CGC shell script's embedded authentication has been replaced with COSMIC_USERNAME/COSMIC_PASSWORD variables; its legacy network workflow is not part of the current validated pipeline and has not been executed.

The bat COSMIC and whole-genome experimental APOBEC calculations have completed on the user's HPC according to the supplied full audit logs. The experimental values still need table export and historical numerical comparison; archive validation alone does not establish agreement with the original biological results.
