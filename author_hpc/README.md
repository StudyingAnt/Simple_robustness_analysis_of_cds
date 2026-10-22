# Author HPC execution archive

The authors ran the full human, bat and experimental-profile calculations and downstream analyses on **Stony Brook University's SeaWulf HPC system**. SeaWulf requires institutional access; this repository does not provide access to that system. Reader instructions are the direct terminal/Python commands in the [reproduction guide](../MANUSCRIPT_REPRODUCTION.md) and each stage's README.

These folders preserve the actual Slurm scripts, recorded paths/job settings, completion checks and pilot material used by the authors. Researchers using an HPC system may consult these examples after adapting scheduler, queue/partition, account, module/Conda setup, CPU/memory/time requests, mail address and filesystem paths to their own site. They are not generic submission commands or prerequisites for the portable workflow.

| Stage | Author execution record |
| --- | --- |
| Human/COSMIC input preparation | [01 author_hpc](../01_data_downloading_and_preprocessing/author_hpc/README.md) |
| Human RUMS, audit and export | [02 author_hpc](../02_run_sinabro/author_hpc/README.md) |
| Signature distributions | [03 author_hpc](../03_mean_std_analysis/author_hpc/README.md) |
| PCA and manuscript analyses | [07 author_hpc](../07_manuscript_validation/author_hpc/README.md) (shared downstream job) |
| Bat extraction and COSMIC | [06 author_hpc](../06_Bat_genome/author_hpc/README.md) |
| Experimental profiles | [experimental author_hpc](../06_Bat_genome/experimental_apobec/author_hpc/README.md) |

## Recorded environment activation

Where required by the site's installed module configuration:

```bash
module load anaconda/3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ABSOLUTE/PATH/TO/YOUR/ENVIRONMENT
export SINABRO_ENV="$CONDA_PREFIX"
export SINABRO_PYTHON="$CONDA_PREFIX/bin/python"
export SINABRO_SOURCE="$(cd ../sinabro && pwd)"
```

Export the environment prefix before submitting the Slurm wrappers, particularly when using a custom installation path. Reference, package and output directories must be visible from compute nodes. The wrappers' standard project layout is described in [the reproduction guide](../MANUSCRIPT_REPRODUCTION.md).


## Layout change for the authors' existing checkout

Slurm files now live at `STAGE/author_hpc/FILE_seawulf.slurm`. Pilot controllers moved to the same author folder; their import paths and wrapper calls are adjusted. Main calculation and export Python entry points retain their original locations and scientific source bytes. Existing data, checkpoints and environment prefixes do not need recreation. Use the new wrapper paths for any future submissions; already submitted/completed jobs do not need rerunning.

Full checkpoint archives and raw/processed restricted COSMIC tables remain external to Git. The checked record describes the authors' reports, not an independent download of their full HPC filesystem. Before this reorganization, the entire checkout was preserved on [backup/before-reader-workflow-2026-10-22](https://github.com/StudyingAnt/Simple_robustness_analysis_of_cds/tree/backup/before-reader-workflow-2026-10-22).
