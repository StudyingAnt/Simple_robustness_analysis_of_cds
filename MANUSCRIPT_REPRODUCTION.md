# Manuscript reproduction map

This repository contains manuscript-specific analyses; the [Sinabro package](https://github.com/StudyingAnt/sinabro) contains the simulator. The corrected implementation and newly generated results are distinguished from historical files. Exact historical random streams and every final figure input were not recovered.

## Execution order

1. Create the [pinned Conda environment](CONDA_SETUP.md) and install the separate package checkout.
2. [Download and preprocess GENCODE v40](01_data_downloading_and_preprocessing/README.md), then [simulate, audit and export human results](02_run_sinabro/README.md).
3. [Prepare the frozen bat cohorts](06_Bat_genome/HPC_REPRODUCTION.md), then [simulate, audit and export bat COSMIC controls](06_Bat_genome/COSMIC_PRODUCTION.md).
4. [Prepare the six experimental APOBEC profiles and run all three species](06_Bat_genome/experimental_apobec/README.md), then audit and export.
5. [Generate the manuscript analyses and independently validate the signature sampler](07_manuscript_validation/README.md). These two jobs are independent once the full exports exist.
6. Run the separately labelled S11/S12 reconstructions where their documented scope is sufficient. Retain the unavailable-input statements when reporting them.

For an existing completed campaign, start from its audited exports rather than repeating steps 2–4. Use the checksums and scientific configuration in the manifests to verify compatibility.

## Figure and method coverage

| Original element | Current entry point / data | Evidence and remaining limits |
|---|---|---|
| Figures 1 and S10: algorithm and conceptual schematics | Separate Sinabro package; methods | Conceptual artwork is not a simulation output. Package tests check computational behavior. |
| Figure 2 and signature rank/distribution panels | `03_mean_std_analysis/analyze_signature_distributions.py` | Full human analysis reported passed in HPC job 2308077; SBS2 leads mean and SD in all three signature groups. New plots use all observations. Historical supplementary artwork is not claimed to be pixel-identical. |
| Figures 3–5 and related PCA, GC3 and signature-feature analyses | `07_manuscript_validation/validate_manuscript_results.py` | Local numerical and bounded end-to-end checks pass. Full new HPC analysis review is pending. An independent historical-table crosscheck recovered the reported rounded PCA variance and transition regression. |
| Figure 6: experimental motif preferences | Six-enzyme source workbook and `06_Bat_genome/experimental_apobec/prepare_profiles.py` | Workbook identities, values and profile conversion checked. Upstream raw assay reads and Galaxy processing have not been rerun; profile preparation is not an independent replication of the wet-lab assay. |
| Figure 7 and S8/S9: species and experimental-enzyme comparisons | Audited human, bat COSMIC and experimental exports; combined downstream analysis | Full trajectory campaigns completed and audited from supplied HPC logs. The new native-process means do not preserve every historical mean ordering; actual statistics and panels must be reviewed. Transcript-level comparisons are not orthologue-matched biological replicates. |
| S1/S2: random-sequence signature validation | `07_manuscript_validation/validate_signature_sampler.py` | Actual sampler checked locally against an independent conditional oracle. Full 79-signature validation review pending. |
| S11: relaxed BLOSUM criterion | [Recovered historical sources](07_manuscript_validation/HISTORICAL_SUPPLEMENTARY_INPUTS.md); separate bounded reconstruction runner | 3,740 historical BLOSUM values and matching transcript sequences recovered. Original SVG compares 73,214 identity points with 3,740 BLOSUM points; a same-transcript comparison is now possible. The twelve-case example is a different, explicitly new subset. |
| S12: BRCA1 graph | Same runner, `--stage s12` | Bounded BRCA1 simulation passes. The reconstruction uses a documented radial layout. Observed COSMIC breast-cancer variant overlay input is missing, so the original red-point panel remains unreproduced. |
| Static supporting tables | `00_supporting_information/` | Original workbooks retained. They must not be relabelled as numerical outputs of a new campaign. Regenerated PCA tables are exported separately. |

The full trajectory campaigns comprise 5,783,906,000 human COSMIC paths, 435,184,000 bat COSMIC paths and 765,672,000 experimental-enzyme paths. Audit success establishes completeness and consistency of saved results; it does not establish every biological conclusion or complete manuscript reproduction.

## Portable code and project layout

Python runners accept explicit paths and do not require SeaWulf. `_seawulf.slurm` wrappers supply scheduler resources, module loading and cluster environment paths. Current downstream analyses run in the pinned Python environment without the historical R environment.

Keep external references in `data/raw/`, prepared inputs in `data/processed/`, parameterized runs/checkpoints and intermediate exports in `runs/`, and final tables/figures in `results/`. Generated data and environment exports are excluded from Git; retain them with the campaign. The six-enzyme workbook is an explicitly included research source.

Preserve each run manifest, full audit report, export manifest, numerical analysis summaries and resolved Conda environment. Report package and analysis commit IDs together. Re-running a seed under changed sampling logic does not establish identity with a historical random stream.

## Release status and historical sources

The current branch is a reproduction candidate. A final paper release must record the remaining full-run validation results and clearly disclose unrecovered supplementary inputs. It is not marked as a completely verified release while those checks are pending.

Original analysis sources are preserved at commit `c2c271e36a0f478fd3b3464902a450c514f5b789` on [the original backup branch](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-cleanup-2026-10-14). Historical scripts and R plotting code retained in this tree provide method provenance; their old fixed paths and package imports are not the current execution instructions. Timing and parallel validation remain research reproducibility checks and are documented separately from the primary run instructions.
