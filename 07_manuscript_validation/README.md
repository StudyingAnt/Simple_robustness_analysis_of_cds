# Manuscript analyses and supplementary panels

These programs use the [pinned Python environment](../CONDA_SETUP.md) and ordinary terminal commands. No HPC account is required. Run from the analysis repository root and provide your own audited upstream exports. The original full calculations were performed on the authors' HPC system; scheduler details are kept separately.

## PCA, sequence/signature features, statistics and species panels

After completing [human table export](../02_run_sinabro/AGGREGATION.md), [bat COSMIC export](../06_Bat_genome/README.md) and [experimental-profile export](../06_Bat_genome/experimental_apobec/README.md):

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
SINABRO_SOURCE=/ABSOLUTE/PATH/TO/YOUR/sinabro

python 07_manuscript_validation/validate_manuscript_results.py \
  --human-export "$PROJECT/runs/human_export" \
  --bat-export "$PROJECT/runs/bat_cosmic_export" \
  --apobec-export "$PROJECT/runs/experimental_apobec_export" \
  --signatures "$SINABRO_SOURCE/sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt" \
  --out-dir "$PROJECT/results/manuscript_analysis" --n-traj 1000
```

Choose a new output directory. The program verifies manifests, populations, signature identity and compatible inputs, then writes `analysis_summary.json`, statistical/feature/PCA tables and PDF/PNG panels. It reuses completed trajectories; it does not rerun the genome simulations. The full manuscript analysis expects the full cohorts, not the synthetic usage example.

PCA uses mean-centered full SVD without standardization and the historical sign convention. All 79-, 60- and 41-signature groups are retained; Figure 4 fits the 60-signature model and displays 41 known-aetiology signatures. Statistical comparisons include Kruskal–Wallis, median-centered Levene, Dunn/Bonferroni and asymptotic Mann–Whitney tests with documented family corrections. Transcript observations are not independent biological replicates. See [methods and reference results](../docs/validation_reference.md).

## Signature sampler figures (S1/S2)

This independent calculation uses random sequences rather than the genome exports:

```bash
python 07_manuscript_validation/validate_signature_sampler.py \
  --sinabro-source "$SINABRO_SOURCE" \
  --out-dir "$PROJECT/results/signature_sampler" \
  --workers 1 --n-sequences 10000 --length 2200 --gc 0.4 --seed 42
```

This samples one event for each of 10,000 random 2,200-bp sequences for all 79 signatures (790,000 events). Outputs include `sampler_validation_summary.json`, channel counts, cosine summaries and figures. The independent conditional probability oracle accounts for available sequence opportunities. A smaller `--n-sequences` or `--signatures` subset is an execution example rather than full S1/S2 reproduction.

## Supplementary analyses

| Analysis | Direct execution guide | Scope |
| --- | --- | --- |
| S11 criterion comparison | [S11_HISTORICAL.md](S11_HISTORICAL.md) | 3,740 matched archived transcript estimates; no new paths |
| S12 BRCA1 graph | [BRCA1_GRAPH.md](BRCA1_GRAPH.md) | Fresh 1,000-path SBS2/BLOSUM80 simulation |
| Reported breast SNVs | [Coordinate preparation](../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md) | COSMIC phenotype selection, genomic mapping and reference checks |
| S12 reported-allele marks | [BRCA1_OBSERVED_SNVS.md](BRCA1_OBSERVED_SNVS.md) | Explicit single-substitution allele matching, using the saved graph |

## Author execution

The authors ran these analyses on Stony Brook University's SeaWulf HPC system, which requires institutional access. [author_hpc/](author_hpc/README.md) preserves the actual Slurm scripts, recorded campaign paths and completion checks. These are examples to adapt to your own site's settings if you use an HPC system; they are not prerequisites for the commands above.
