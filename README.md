# Human and bat genome robustness under COSMIC mutational signatures

This repository contains the manuscript-specific analysis scripts. The separate [Sinabro package](https://github.com/StudyingAnt/sinabro) provides the sequential mutation simulator.

Start with [Conda environment setup](CONDA_SETUP.md), using the pinned [environment.yml](environment.yml). For the full human simulation on SeaWulf, see [Slurm submission and resume](02_run_sinabro/HPC_SUBMISSION.md).

The [manuscript reproduction map](MANUSCRIPT_REPRODUCTION.md) connects the execution order, original figures, numerical checks and remaining source-data limits.

The [portable runner and parallel validation](02_run_sinabro/REPRODUCTION.md) run without Slurm. Python filenames describe the analysis; cluster names are reserved for scheduler-specific wrappers. The separate [small SeaWulf validation job](02_run_sinabro/validate_parallel_seawulf.slurm) compares one and eight workers before full submission.

## Reproduction status

Reproduction is checked one stage at a time. The first stage, **GENCODE v40 download and CDS-plus preprocessing**, has been rerun from a fresh download. The new CDS-plus and fifteen historical rejection files match the preserved original preprocessing script byte for byte.

| Stage | Status |
|---|---|
| GENCODE v40 download and gzip integrity | Verified |
| CDS-plus generation and rejection accounting | Verified |
| Length-stratified pilot: 12 transcripts × 3 signatures × 20 paths | Verified; 720/720 complete |
| Same pilot with 100 paths per combination | Verified; 3,600/3,600 complete; prior 720 paths reproduced |
| Portable one-worker vs eight-worker smoke test | Verified locally; user-reported SeaWulf comparison also passed (3,600 paths per mode) |
| Full human RUMS calculation and saved-output audit | User-reported full audit passed: 73,214 transcripts × 79 signatures × 1,000 paths, no incomplete paths |
| Audited NPZ-to-table export | User-reported full export passed: 73,214 checksums matched; 19,309 genes and 73 internal-stop transcripts retained |
| Descriptive signature ranks and distribution figures | Full HPC analysis passed (job 2308077); SBS2 leads mean and SD in all three signature groups |
| Two-bat COSMIC calculation, audit and export | Passed: 31,185 Molossus and 23,213 Myotis transcripts × eight signatures × 1,000 paths |
| Six experimental enzyme profiles on all three species | Full audit passed (job 2336448), 765,672,000 paths; 18 exported summary combinations checked |
| Original PCA/features/statistical tables and species panels | Full HPC numerical report passed; [results and manuscript implications](07_manuscript_validation/VALIDATION_RESULTS.md). Figure 4's numerical loading comparison is complete; final artwork review remains. |
| Original S1/S2 random-sequence signature validation | Full 79-signature validation passed: 790,000 substitutions, all conditional-oracle checks passed; mean original-profile cosine 0.97519. |
| Original S11 BLOSUM subset and S12 BRCA1 example | [S11 historical matched-population analysis](07_manuscript_validation/S11_HISTORICAL.md) passed on HPC (job 2345768), all 3,740 input sequences verified. S12 observed-variant mapping remains unresolved; fresh sensitivity simulations are separate. |

The complete whole-genome checkpoints and exports are available from the checked human, bat COSMIC and experimental-enzyme runs. Preserve their audit reports and [exported tables](02_run_sinabro/AGGREGATION.md). Continue with [original manuscript downstream and sampler validation](07_manuscript_validation/README.md), which reuses those calculations and reports the actual observed results without requiring the historical conclusions to hold. No claim of full manuscript reproduction is made while the remaining rows above are pending.

Follow [the first-stage commands and checks](01_data_downloading_and_preprocessing/README.md). Reference counts and hashes are stored in [expected_gencode_v40.json](01_data_downloading_and_preprocessing/expected_gencode_v40.json).

The original analysis code is preserved at commit `c2c271e36a0f478fd3b3464902a450c514f5b789` on [backup/before-cleanup-2026-10-14](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-cleanup-2026-10-14). Later-stage scripts still reflect the historical workflow and may require path and implementation corrections before execution. Preprocessing agreement alone does not establish reproduction of the manuscript figures.

Downloaded input and generated output files under `data/` are excluded from Git. Each first-stage run generates manifests with counts and hashes. Historical results and outputs produced using a corrected Sinabro implementation must be distinguished during numerical validation.

BSD-3-Clause; see [LICENSE](LICENSE).
