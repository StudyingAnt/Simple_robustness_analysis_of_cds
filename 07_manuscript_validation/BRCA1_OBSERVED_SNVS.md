# BRCA1 graph annotation with reported breast-cancer SNVs

The author-reported COSMIC preparation has passed: **354 breast BRCA1 observations**, of which **235** map as SNVs to the simulated CDS; these represent **174 distinct target substitutions**. There are 74 non-SNV/unknown-allele exclusions and 45 outside-target-CDS exclusions. The input archive, FASTA, map and processing-code hashes were verified against the local reference files. There were no unresolved phenotypes, reference mismatches or genomic-annotation conflicts. See [coordinate preparation](../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md).

This next job compares the checked SNV set with the existing full graph. It performs **no new simulation, download or package installation**.

## Submit one job on SeaWulf

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/annotate_brca1_graph_seawulf.slurm
```

Allocation: one CPU, 4 GB, 30 minutes. The existing `PROJECT/sinabro-revision` environment is activated.

The Python program discovers successful verified inputs under:

```text
PROJECT/runs/brca1_sbs2_blosum20_n1000_*/
PROJECT/data/processed/cosmic_v101_grch38/brca1_breast_*/
```

It checks the manifests and input-file hashes. Multiple equivalent runs are acceptable if the scientific input files have identical hashes. If scientifically different successful inputs are found, it stops and requires an explicit selection. To choose known directories:

```bash
SINABRO_BRCA1_GRAPH_DIR=/ABSOLUTE/PATH/TO/EXISTING_GRAPH \
COSMIC_BREAST_DIR=/ABSOLUTE/PATH/TO/EXISTING_MAPPED_VARIANTS \
  sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/annotate_brca1_graph_seawulf.slurm
```

No previous successful job needs resubmission. Source directories are read-only inputs to this program. The figure layout and simulation counts remain those of the checked full graph.

## Matching definition

**The red figure marks require exactly one reference-relative CDS substitution in the simulated node, matching a prepared COSMIC breast SNV in position, reference and alternate base.** Matching a position with a different alternate allele is insufficient. The graph reference and every observed variant's target reference base are checked again.

This is allele agreement, not proof that an entire simulated DNA sequence was observed in a tumour. Even a reference-plus-one-SNV state cannot establish the absence of additional unreported tumour variants. All output rows therefore retain `observed_genotype_confirmed=false`.

For nodes with several substitutions, the audit table separately records which component alleles have been reported and whether all component alleles individually occur in the catalog. This does not establish co-occurrence within one sample, and such nodes are not coloured as observed genotypes. The reference node is not called an observed-variant match. Zero exact single-SNV matches is a legitimate outcome and is not replaced with artificial highlights.

The original observation-to-node rule and exact historical COSMIC release remain unrecovered. The regenerated panel uses the explicit rule above and carries a matching caption. It does not establish driver status, causation by APOBEC or experimental functional loss. BLOSUM80 distance remains an alignment-score proxy.

## Outputs and result check

Outputs go to:

```text
PROJECT/results/brca1_cosmic_snv_matches_JOBID/
```

- **`s12_observed_annotation_summary.json`**: input paths/hashes, rerun graph audit, number of matching variants/nodes and annotation status.
- `s12_annotated_nodes.tsv`: graph states with exact-single-SNV/component-match columns and explicit genotype-evidence status.
- `s12_observed_variant_matches.tsv`: every mapped SNV, source observation/sample counts and numbers of graph matches.
- `figureS12_cosmic_single_snv_matches.pdf/.png`: regenerated graph; red points follow the single-SNV allele rule and the black star marks the reference. Other node colours show BLOSUM80 alignment-score distance.
- `figureS12_caption.txt`: standalone caption stating the matching definition and evidence limits.

Share the summary JSON. `annotation_passed=true`, `regenerated_allele_match_panel_created=true`, `bounded_software_check=false` and `fresh_simulation=false` describe successful full-data postprocessing. `original_figure_reproduced=false`, `historical_observation_rule_recovered=false` and the legacy `observed_cosmic_variant_overlay_reproduced=false` remain explicit because the historical panel's matching semantics have not been recovered. These are distinct from successful generation of the newly defined allele-match panel.

The summary retains `observed_genotypes_confirmed=false`. The completed numerical/simulation analyses are unaffected. The real-data overlap counts must be read from this new HPC result before making a statement about how many of the 174 variants match. This is separate from the local plotting check, which used three synthetic alleles rather than actual COSMIC observations.

## Validation and portable execution

Local validation rechecked the actual full graph: 1,000 complete paths, 3,152 reconstructed sequence states, 3,623 single-substitution edges. Ten dedicated tests cover exact/component semantics, alternate-allele mismatches, manifest integrity, source selection, failed preparation, zero-match plotting and output preservation. Only synthetic observations were used in these annotation tests; restricted COSMIC input tables remain on HPC.

An end-to-end local software check also generated annotated tables, summary, caption and images from that full graph with three explicitly synthetic alleles. Its report is labelled `bounded_software_check=true`, and its plot/caption say SOFTWARE CHECK. Its three matches are not a result for the actual 174 COSMIC variants.

On another computer with the pinned analysis environment, run from the analysis checkout:

```bash
python 07_manuscript_validation/annotate_brca1_graph.py \
  --graph-dir /PATH/TO/EXISTING_FULL_GRAPH \
  --variants-dir /PATH/TO/EXISTING_PREPARED_VARIANTS \
  --out-dir /PATH/TO/NEW_RESULTS_DIRECTORY
```

Slurm settings remain isolated in the SeaWulf wrapper. The annotation uses the existing graph audit code and NumPy/Biopython/Matplotlib/SciPy environment, without running Sinabro trajectories.
