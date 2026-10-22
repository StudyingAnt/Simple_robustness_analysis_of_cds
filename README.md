# Human and bat genome robustness under COSMIC mutational signatures

Analysis code and reproducible workflows accompanying *Human and bats genome robustness under COSMIC mutational signatures*. The separate [Sinabro package](https://github.com/StudyingAnt/sinabro) provides the sequential mutation simulator; this repository provides reference preparation, human/bat calculations, statistical analyses and figure generation.

## Start here

1. [Install the pinned environment and Sinabro](CONDA_SETUP.md).
2. Follow the [manuscript reproduction guide](MANUSCRIPT_REPRODUCTION.md) in calculation order.
3. Consult [data sources and availability](DATA_AVAILABILITY.md), [methods and validation](docs/validation_reference.md), and the [software manifest](software_manifest.json).

Python workflows run on other computers using explicit paths. Files ending in `_seawulf.slurm` contain SeaWulf scheduler settings. Full calculations can be distributed over job-array shards and worker processes; smaller examples and integrity checks are provided before full submission.

## Analysis workflows

| Stage | Entry point |
| --- | --- |
| Human reference download and CDS-plus preparation | [GENCODE v40](01_data_downloading_and_preprocessing/README.md) |
| Human RUMS, checkpoint audit and table export | [Human calculation](02_run_sinabro/README.md) |
| Signature distributions and summary ranks | [Distribution analysis](03_mean_std_analysis/README.md) |
| PCA, sequence features and manuscript statistics | [Manuscript analyses](07_manuscript_validation/README.md) |
| Bat CDS cohorts and COSMIC controls | [Bat workflows](06_Bat_genome/README.md) |
| Six experimental APOBEC profiles across three species | [Experimental profiles](06_Bat_genome/experimental_apobec/README.md) |
| Supplementary BLOSUM criterion comparison | [Archived S11 matched-transcript analysis](07_manuscript_validation/S11_HISTORICAL.md) |
| BRCA1 mutation graph and reported SNV matches | [S12 graph](07_manuscript_validation/BRCA1_GRAPH.md) and [allele annotation](07_manuscript_validation/BRCA1_OBSERVED_SNVS.md) |

## Reproducibility

The checked workflow uses **Python 3.9.15, NumPy 1.21.6, Biopython 1.79, SciPy 1.8.1 and Matplotlib 3.6.2**, with **Sinabro 3.0.0**. The reference campaigns comprise 5,783,906,000 human COSMIC paths, 435,184,000 bat COSMIC paths and 765,672,000 experimental-profile paths. Saved-output audits reported no incomplete paths. The [validation reference](docs/validation_reference.md) documents numerical comparisons, sampler checks, source fingerprints and the scope of supplementary reconstructions.

The present edition contains corrected package workflows. Historical random streams are not recovered; archived S11 estimates and the newly defined S12 allele-match panel are labelled explicitly. The [historical source archive](historical/README.md) preserves earlier scripts and methodological provenance separately from the current execution instructions.

## Project data layout

Keep generated files outside the source checkout:

```text
PROJECT/
  sinabro/                              # package checkout
  Simple_robustness_analysis_of_cds/     # this checkout
  data/raw/                             # external source files
  data/processed/                       # prepared inputs and variant tables
  runs/                                 # parameterized calculations/checkpoints
  results/                              # final tables and figures
  logs/                                 # scheduler output
```

Raw references, COSMIC observations, generated outputs and environment exports are not distributed in Git. Frozen cohort IDs, historical S11 estimates, the coordinate reference and the experimental preference workbook are included where required; see [data availability](DATA_AVAILABILITY.md).

## Tests and citation

With the environment and package installed:

```bash
python -m unittest discover -s tests -v
```

Tests cover preparation, deterministic sampling/parallelism, checkpoint audits, exports, statistics and supplementary matching rules. Large reference campaigns are separate from software tests. [CITATION.cff](CITATION.cff) provides software attribution; cite the manuscript and record both repository commits with scientific outputs.

BSD-3-Clause. See [LICENSE](LICENSE).
