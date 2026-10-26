# Manuscript reproduction guide

This repository supplies analysis workflows for *Human and bats genome robustness under COSMIC mutational signatures*. [Sinabro](https://github.com/StudyingAnt/sinabro) supplies the simulator. Start with [environment setup](CONDA_SETUP.md), [data availability](DATA_AVAILABILITY.md) and the [software manifest](software_manifest.json).

## Execution order

| Step | Input and calculation | Output |
| --- | --- | --- |
| 1 | [GENCODE v40 download and CDS-plus preparation](01_data_downloading_and_preprocessing/README.md) | 73,214 human CDS-plus records and provenance manifest |
| 2 | [Human COSMIC calculation](02_run_sinabro/README.md), [audit and export](02_run_sinabro/AGGREGATION.md) | Checkpoints, per-transcript results and signature summaries |
| 3 | [Signature distributions](03_mean_std_analysis/README.md) | Rank tables, full-data distributions and plots |
| 4 | [Bat reference preparation and COSMIC controls](06_Bat_genome/README.md) | Frozen Molossus/Myotis cohorts and audited results |
| 5 | [Six experimental APOBEC profiles](06_Bat_genome/experimental_apobec/README.md) | Profile tables and audited results on all three species |
| 6 | [Manuscript statistics and sampler checks](07_robustness_analysis/README.md) | PCA, GC3/signature features, statistical tables and panels |
| 7 | [Historical S11 matched-transcript comparison](07_robustness_analysis/S11_HISTORICAL.md) | Postprocessing of 3,740 archived paired estimates |
| 8 | [BRCA1 graph](07_robustness_analysis/BRCA1_GRAPH.md), [COSMIC breast-SNV preparation](01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md), [allele annotation](07_robustness_analysis/BRCA1_OBSERVED_SNVS.md) | Audited 1,000-path graph, observed-allele matches and supplementary panel |

Run each linked stage from the analysis repository root using its direct Python/terminal commands. To try the workflow without genome downloads, start with the [included synthetic example](02_run_sinabro/README.md#try-the-workflow-on-the-included-example). Full manuscript campaigns are large calculations. Existing audited exports can be reused to regenerate plots.

## Figure and method map

| Manuscript element | Current workflow | Scope |
| --- | --- | --- |
| Figures 1 / S10 | Package algorithm and conceptual artwork | Software behavior is checked by package tests |
| Figure 2 and rank/distribution panels | `03_mean_std_analysis/analyze_signature_distributions.py` | All human transcript observations and preserved signature groups |
| Figures 3–5 and related panels | `07_robustness_analysis/analyze_robustness_results.py` | Mean-centered PCA, GC3, signature features and original statistical comparisons |
| Figure 6 | Six-enzyme workbook and `prepare_profiles.py` | Motif-average preferences; upstream wet-lab/read processing is outside this workflow |
| Figure 7 / S8 / S9 | Human, bat and experimental-profile exports | Transcript-level comparisons; not independent biological replicates or an orthologue-matched design |
| S1 / S2 | `validate_signature_sampler.py` | 79 signatures, 790,000 single-substitution samples and an independent conditional oracle |
| S11 | `analyze_s11_historical.py` | Matched archived estimates, not a fresh corrected-package simulation of the full cohort |
| S12 | `run_supplementary_cases.py`, coordinate preparation and `annotate_brca1_graph.py` | Fresh simulations with an explicit current single-SNV allele-matching rule |
| Supporting tables | `00_supporting_information/` | Preserved source workbooks; newly generated tables are separate outputs |

See [methods and reference results](docs/validation_reference.md) for quantitative agreement and differences. The current computations use corrected package behavior and explicit seeds. Exact historical random streams are unavailable, and historical estimates are not relabelled as fresh simulation outputs.

## Your project directory

```text
PROJECT/
  sinabro/
  Simple_robustness_analysis_of_cds/
  data/raw/
  data/processed/
  runs/
  results/
  logs/
```

Use an absolute `PROJECT` path on your computer and keep generated data outside the source checkout. Each guide names the inputs and expected output files. Choose the number of local processes with `--workers`; one process is valid. Preserve manifests, input hashes, resolved environments and both repository commits with each calculation.

The authors used Stony Brook University's SeaWulf HPC system for full calculations. It requires institutional access and is not part of the reader setup. Each stage's `author_hpc/` folder preserves actual submitted scripts, reference campaign paths, completion records and pilot materials. If you use an HPC system, adapt those examples to your own site's settings. See the [author execution index](author_hpc/README.md); these records are optional for running the portable workflows.

## Supplementary interpretation and historical sources

For S11, the recovered 3,740 estimates are paired to the same transcript population under identity and BLOSUM criteria. The original plotted populations were unequal; the current postprocessing uses matched populations and retains archived provenance.

For S12, genomic coordinates are mapped to BRCA1-203, with strand and reference-base checks. Red nodes denote single-substitution states matching reported COSMIC v101 breast SNVs. Individually reported component alleles do not establish complete tumour genotypes. The original source release/node-matching rule is not recovered; the regenerated caption states the current definition.

Earlier Python/R/SGE sources are in [historical/](historical/README.md); author execution and validation records are indexed in [author_hpc/](author_hpc/README.md), separate from the current workflow. Original package implementations and the pre-publication analysis snapshot remain in backup branches. This repository is code for manuscript reproduction; journal publication metadata and a persistent archival DOI should be added when available rather than inferred from a repository name.
