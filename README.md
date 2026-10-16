# Human and bat genome robustness under COSMIC mutational signatures

This repository contains the manuscript-specific analysis scripts. The separate [Sinabro package](https://github.com/StudyingAnt/sinabro) provides the sequential mutation simulator.

Start with [Conda environment setup](CONDA_SETUP.md), using the pinned [environment.yml](environment.yml). For the full human simulation on SeaWulf, see [Slurm submission and resume](02_run_sinabro/HPC_SUBMISSION.md).

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
| Audited NPZ-to-table export | Implemented and checked on a small real local run; full HPC export pending |
| Downstream statistics, PCA and figures | Pending reproduction |
| Bat data and analyses | Pending validation |

The full human calculation and archive audit were reported complete on 2026-10-16 (audit job 2304985). Next, use [the audited table export](02_run_sinabro/AGGREGATION.md), preserving the checkpoint files and audit reports. The export prepares inputs for the historical analysis scripts; reproducing manuscript conclusions and figures remains a separate stage. See [pilot commands and instrumentation](02_run_sinabro/README.md) for the earlier bounded measurements.

Follow [the first-stage commands and checks](01_data_downloading_and_preprocessing/README.md). Reference counts and hashes are stored in [expected_gencode_v40.json](01_data_downloading_and_preprocessing/expected_gencode_v40.json).

The original analysis code is preserved at commit `c2c271e36a0f478fd3b3464902a450c514f5b789` on [backup/before-cleanup-2026-10-14](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-cleanup-2026-10-14). Later-stage scripts still reflect the historical workflow and may require path and implementation corrections before execution. Preprocessing agreement alone does not establish reproduction of the manuscript figures.

Downloaded input and generated output files under `data/` are excluded from Git. Each first-stage run generates manifests with counts and hashes. Historical results and outputs produced using a corrected Sinabro implementation must be distinguished during numerical validation.

BSD-3-Clause; see [LICENSE](LICENSE).
