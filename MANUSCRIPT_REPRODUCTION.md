# Manuscript reproduction guide

This repository supplies analysis workflows for *Human and bats genome robustness under COSMIC mutational signatures*. [Sinabro](https://github.com/StudyingAnt/sinabro) supplies the simulator. Start with [environment setup](CONDA_SETUP.md), [data availability](DATA_AVAILABILITY.md) and the [software manifest](software_manifest.json).

## Execution order

| Step | Input and calculation | Output |
| --- | --- | --- |
| 1 | [GENCODE v40 download and CDS-plus preparation](01_data_downloading_and_preprocessing/README.md) | 73,214 human CDS-plus records and provenance manifest |
| 2 | [Human COSMIC calculation](02_run_sinabro/README.md), [Slurm submission/resume](02_run_sinabro/HPC_SUBMISSION.md), [audit and export](02_run_sinabro/AGGREGATION.md) | Checkpoints, per-transcript results and signature summaries |
| 3 | [Signature distributions](03_mean_std_analysis/README.md) | Rank tables, full-data distributions and plots |
| 4 | [Bat reference preparation](06_Bat_genome/HPC_REPRODUCTION.md) and [COSMIC controls](06_Bat_genome/COSMIC_PRODUCTION.md) | Frozen Molossus/Myotis cohorts and audited results |
| 5 | [Six experimental APOBEC profiles](06_Bat_genome/experimental_apobec/README.md) | Profile tables and audited results on all three species |
| 6 | [Manuscript statistics and sampler checks](07_manuscript_validation/README.md) | PCA, GC3/signature features, statistical tables and panels |
| 7 | [Historical S11 matched-transcript comparison](07_manuscript_validation/S11_HISTORICAL.md) | Postprocessing of 3,740 archived paired estimates |
| 8 | [BRCA1 graph](07_manuscript_validation/BRCA1_GRAPH.md), [COSMIC breast-SNV preparation](01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md), [allele annotation](07_manuscript_validation/BRCA1_OBSERVED_SNVS.md) | Audited 1,000-path graph, observed-allele matches and supplementary panel |

Start new calculations with the small parallel check in [REPRODUCTION.md](02_run_sinabro/REPRODUCTION.md). Existing audited campaigns can proceed from their exported tables. Do not resimulate completed trajectories merely to regenerate plots.

## Figure and method map

| Manuscript element | Current workflow | Scope |
| --- | --- | --- |
| Figures 1 / S10 | Package algorithm and conceptual artwork | Software behavior is checked by package tests |
| Figure 2 and rank/distribution panels | `03_mean_std_analysis/analyze_signature_distributions.py` | All human transcript observations and preserved signature groups |
| Figures 3–5 and related panels | `07_manuscript_validation/validate_manuscript_results.py` | Mean-centered PCA, GC3, signature features and original statistical comparisons |
| Figure 6 | Six-enzyme workbook and `prepare_profiles.py` | Motif-average preferences; upstream wet-lab/read processing is outside this workflow |
| Figure 7 / S8 / S9 | Human, bat and experimental-profile exports | Transcript-level comparisons; not independent biological replicates or an orthologue-matched design |
| S1 / S2 | `validate_signature_sampler.py` | 79 signatures, 790,000 single-substitution samples and an independent conditional oracle |
| S11 | `analyze_s11_historical.py` | Matched archived estimates, not a fresh corrected-package simulation of the full cohort |
| S12 | `run_supplementary_cases.py`, coordinate preparation and `annotate_brca1_graph.py` | Fresh simulations with an explicit current single-SNV allele-matching rule |
| Supporting tables | `00_supporting_information/` | Preserved source workbooks; newly generated tables are separate outputs |

See [methods and reference results](docs/validation_reference.md) for quantitative agreement and differences. The current computations use corrected package behavior and explicit seeds. Exact historical random streams are unavailable, and historical estimates are not relabelled as fresh simulation outputs.

## Project layout and scheduler use

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

Submit Slurm files from `PROJECT/`, with `logs/` created first. `PROJECT`, `ANALYSIS_DIR`, `SINABRO_ENV`, `SINABRO_SOURCE`, input and output variables are configurable in the wrappers. Current executable paths are retained so existing HPC workflows remain valid. Python programs accept explicit paths without requiring Slurm.

The completed example campaign locations in the wrappers identify the reference calculations. On another project, override the corresponding variables/arguments to point at that project's exports. Raw/processed data and generated figures belong outside Git; preserve each manifest, complete audit, input hash, resolved environment and both repository commits with the campaign.

## Supplementary interpretation and historical sources

For S11, the recovered 3,740 estimates are paired to the same transcript population under identity and BLOSUM criteria. The original plotted populations were unequal; the current postprocessing uses matched populations and retains archived provenance.

For S12, genomic coordinates are mapped to BRCA1-203, with strand and reference-base checks. Red nodes denote single-substitution states matching reported COSMIC v101 breast SNVs. Individually reported component alleles do not establish complete tumour genotypes. The original source release/node-matching rule is not recovered; the regenerated caption states the current definition.

Earlier Python/R/SGE sources and detailed validation records are in [historical/](historical/README.md), separate from the current workflow. Original package implementations and the pre-publication analysis snapshot remain in backup branches. This repository is code for manuscript reproduction; journal publication metadata and a persistent archival DOI should be added when available rather than inferred from a repository name.
