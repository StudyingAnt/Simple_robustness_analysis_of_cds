# Human reference and observed-variant preparation: author HPC record

The authors ran this stage on Stony Brook University's **SeaWulf** HPC system, which requires institutional access. For execution on your own computer, start with the [stage README](../README.md). These are the authors' execution examples, not a requirement for readers. If you use your own HPC system, adapt the scheduler, site/environment settings and paths before use.

The initial GENCODE download/preprocessing jobs were submitted using preambles shared interactively; standalone copies of those initial Slurm files are not present in this repository. The portable scripts and recorded input fingerprints are retained. No replacement script is presented as an original submitted file.

## Execution guides and records

- [Actual COSMIC retrieval/source-extraction settings and full-data report](COSMIC_VARIANTS.md)
- [Actual breast-SNV preparation and coordinate-mapping record](COSMIC_BREAST_VARIANTS.md)

## Retained Slurm files

- [prepare_cosmic_brca1_seawulf.slurm](prepare_cosmic_brca1_seawulf.slurm)
- [prepare_cosmic_breast_variants_seawulf.slurm](prepare_cosmic_breast_variants_seawulf.slurm)

Commands in the execution records are issued from the authors' project root. Full raw/checkpoint data are external to Git. See the [author archive index](../../author_hpc/README.md) for the site/environment context and layout change.
