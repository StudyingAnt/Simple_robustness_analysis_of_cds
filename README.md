# Human and bat genome robustness under COSMIC mutational signatures

This repository contains the manuscript-specific analysis scripts. The separate [Sinabro package](https://github.com/StudyingAnt/sinabro) provides the sequential mutation simulator.

## Reproduction status

Reproduction is checked one stage at a time. The first stage, **GENCODE v40 download and CDS-plus preprocessing**, has been rerun from a fresh download. The new CDS-plus and fifteen historical rejection files match the preserved original preprocessing script byte for byte.

| Stage | Status |
|---|---|
| GENCODE v40 download and gzip integrity | Verified |
| CDS-plus generation and rejection accounting | Verified |
| Length-stratified pilot: 12 transcripts × 3 signatures × 20 paths | Verified; 720/720 complete |
| Full human RUMS calculation, statistics, figures | Pending validation |
| Bat data and analyses | Pending validation |

The next decision is whether to expand the pilot to 100 trajectories per combination. No such run has been performed. See [pilot commands and instrumentation](02_run_sinabro/README.md).

Follow [the first-stage commands and checks](01_data_downloading_and_preprocessing/README.md). Reference counts and hashes are stored in [expected_gencode_v40.json](01_data_downloading_and_preprocessing/expected_gencode_v40.json).

The original analysis code is preserved at commit `c2c271e36a0f478fd3b3464902a450c514f5b789` on [backup/before-cleanup-2026-10-14](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-cleanup-2026-10-14). Later-stage scripts still reflect the historical workflow and may require path and implementation corrections before execution. Preprocessing agreement alone does not establish reproduction of the manuscript figures.

Downloaded input and generated output files under `data/` are excluded from Git. Each first-stage run generates manifests with counts and hashes. Historical results and outputs produced using a corrected Sinabro implementation must be distinguished during numerical validation.

BSD-3-Clause; see [LICENSE](LICENSE).
