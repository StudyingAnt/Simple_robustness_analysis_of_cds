# S12: annotate reported breast-SNV alleles

First complete the [BRCA1 graph](BRCA1_GRAPH.md) and [COSMIC breast-SNV mapping](../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md). This step reads those saved outputs and runs no new trajectories. From the analysis repository root:

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
python 07_manuscript_validation/annotate_brca1_graph.py \
  --graph-dir "$PROJECT/runs/brca1_sbs2_blosum20_n1000" \
  --variants-dir "$PROJECT/data/processed/cosmic_v101_grch38/brca1_breast" \
  --out-dir "$PROJECT/results/brca1_cosmic_snv_matches"
```

Choose a new output directory. Read `s12_observed_annotation_summary.json`; successful processing reports `annotation_passed=true` and `regenerated_allele_match_panel_created=true`. Outputs include annotated states, variant matches, PDF/PNG panel and `figureS12_caption.txt`.

## Matching rule

Red marks require **exactly one reference-relative CDS substitution**, matching a prepared reported SNV in position, reference base and alternate base. Matching a position alone is insufficient. Multi-substitution nodes with reported component alleles are recorded separately and are not red-marked as observed tumour genotypes. Individual catalogued alleles do not establish complete sequence/genotype co-occurrence; `observed_genotypes_confirmed=false` is retained.

Reference full-data results: 174 distinct mapped breast SNVs (235 source observations), **19** matching single-substitution nodes, **288** multi-substitution nodes with a reported component, and 307 nodes with any reported component. These overlap counts are descriptive, not enrichment, mutation-probability or tumour-prevalence estimates.

The historical observation-to-node rule and exact original COSMIC release remain unresolved. Accordingly `original_figure_reproduced=false` and the legacy `observed_cosmic_variant_overlay_reproduced=false` distinguish the newly defined allele-match panel from exact historical reproduction. They do not negate successful annotation under the declared rule.

## Caption definition

Simulated BRCA1-203 (ENST00000357654.9) paths under SBS2. Nodes represent DNA states from 1,000 trajectories stopped at a BLOSUM80 alignment-score distance of at least 20. Colour shows BLOSUM80 distance and radius shows reference DNA Hamming distance. Red nodes have exactly one reference-relative substitution matching a COSMIC v101 SNV with primary site breast (n=19). These allele matches do not establish complete tumour genotypes.

The authors performed this stage on SeaWulf. Their [actual submission and verified output record](author_hpc/BRCA1_OBSERVED_SNVS.md) are retained separately. SeaWulf requires institutional access; the direct command above does not. Adapt the optional HPC example to your own site as needed.
