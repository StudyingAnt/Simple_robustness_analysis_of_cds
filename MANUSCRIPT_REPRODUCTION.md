# Manuscript reproduction map

This repository contains manuscript-specific analyses; the [Sinabro package](https://github.com/StudyingAnt/sinabro) contains the simulator. The corrected implementation and newly generated results are distinguished from historical files. Exact historical random streams and every final figure input were not recovered.

## Execution order

1. Create the [pinned Conda environment](CONDA_SETUP.md) and install the separate package checkout.
2. [Download and preprocess GENCODE v40](01_data_downloading_and_preprocessing/README.md), then [simulate, audit and export human results](02_run_sinabro/README.md).
3. [Prepare the frozen bat cohorts](06_Bat_genome/HPC_REPRODUCTION.md), then [simulate, audit and export bat COSMIC controls](06_Bat_genome/COSMIC_PRODUCTION.md).
4. [Prepare the six experimental APOBEC profiles and run all three species](06_Bat_genome/experimental_apobec/README.md), then audit and export.
5. [Generate the manuscript analyses and independently validate the signature sampler](07_manuscript_validation/README.md). These two jobs are independent once the full exports exist.
6. Run the [historical S11 matched-transcript analysis](07_manuscript_validation/S11_HISTORICAL.md) with its bundled archived estimates and pinned human input. Treat new S11/S12 sensitivity simulations separately; retain the S12 unavailable-observed-input statements.

For an existing completed campaign, start from its audited exports rather than repeating steps 2–4. Use the checksums and scientific configuration in the manifests to verify compatibility.

## Figure and method coverage

| Original element | Current entry point / data | Evidence and remaining limits |
|---|---|---|
| Figures 1 and S10: algorithm and conceptual schematics | Separate Sinabro package; methods | Conceptual artwork is not a simulation output. Package tests check computational behavior. |
| Figure 2 and signature rank/distribution panels | `03_mean_std_analysis/analyze_signature_distributions.py` | Full human analysis reported passed in HPC job 2308077; SBS2 leads mean and SD in all three signature groups. New plots use all observations. Historical supplementary artwork is not claimed to be pixel-identical. |
| Figures 3–5 and related PCA, GC3 and signature-feature analyses | `07_manuscript_validation/validate_manuscript_results.py` | Local numerical and bounded end-to-end checks pass. Full HPC numerical analysis passed and preserves rounded PCA variance and transition-regression values. Figure 4 numerical loading comparison is complete; final artwork review remains; see [validation results](07_manuscript_validation/VALIDATION_RESULTS.md). |
| Figure 6: experimental motif preferences | Six-enzyme source workbook and `06_Bat_genome/experimental_apobec/prepare_profiles.py` | Workbook identities, values and profile conversion checked. Upstream raw assay reads and Galaxy processing have not been rerun; profile preparation is not an independent replication of the wet-lab assay. |
| Figure 7 and S8/S9: species and experimental-enzyme comparisons | Audited human, bat COSMIC and experimental exports; combined downstream analysis | Full trajectory campaigns completed and audited from supplied HPC logs. Full downstream numerical report passed. Myotis has the highest native-panel median and the original rank-test directions are retained, but the human mean is slightly higher; describe the statistic explicitly. Final artwork review remains. Transcript-level comparisons are not orthologue-matched biological replicates. |
| S1/S2: random-sequence signature validation | `07_manuscript_validation/validate_signature_sampler.py` | Full 79-signature / 790,000-substitution validation passed all conditional-oracle checks. Mean original-profile cosine 0.97519; see [validation results](07_manuscript_validation/VALIDATION_RESULTS.md). |
| S11: relaxed BLOSUM criterion | [Historical matched-transcript analysis](07_manuscript_validation/S11_HISTORICAL.md) | All 3,740 historical pairs verified locally against the pinned sequences. Portable code and bundled estimates reproduce matched-population postprocessing. Fresh corrected-package BLOSUM simulation is separate; the twelve-case example is not the original cohort. |
| S12: BRCA1 graph | Same runner, `--stage s12` | Bounded BRCA1 simulation passes. The reconstruction uses a documented radial layout. Observed COSMIC breast-cancer variant overlay input is missing, so the original red-point panel remains unreproduced. |
| Static supporting tables | `00_supporting_information/` | Original workbooks retained. They must not be relabelled as numerical outputs of a new campaign. Regenerated PCA tables are exported separately. |

The full trajectory campaigns comprise 5,783,906,000 human COSMIC paths, 435,184,000 bat COSMIC paths and 765,672,000 experimental-enzyme paths. Audit success establishes completeness and consistency of saved results; it does not establish every biological conclusion or complete manuscript reproduction.

## Portable code and project layout

Python runners accept explicit paths and do not require SeaWulf. `_seawulf.slurm` wrappers supply scheduler resources, module loading and cluster environment paths. Current downstream analyses run in the pinned Python environment without the historical R environment.

Keep external references in `data/raw/`, prepared inputs in `data/processed/`, parameterized runs/checkpoints and intermediate exports in `runs/`, and final tables/figures in `results/`. Generated data and environment exports are excluded from Git; retain them with the campaign. The six-enzyme workbook is an explicitly included research source.

Preserve each run manifest, full audit report, export manifest, numerical analysis summaries and resolved Conda environment. Report package and analysis commit IDs together. Re-running a seed under changed sampling logic does not establish identity with a historical random stream.

## Release status and historical sources

The current branch is a reproduction candidate. Both core full-run numerical reports have passed. A final paper release must finish artwork review, record the supplementary scope and clearly disclose unrecovered observed-variant inputs. It is not marked as a completely verified release while those items remain.

Original analysis sources are preserved at commit `c2c271e36a0f478fd3b3464902a450c514f5b789` on [the original backup branch](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-cleanup-2026-10-14). Historical scripts and R plotting code retained in this tree provide method provenance; their old fixed paths and package imports are not the current execution instructions. Timing and parallel validation remain research reproducibility checks and are documented separately from the primary run instructions.
