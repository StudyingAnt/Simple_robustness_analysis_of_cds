# Manuscript analyses and supplementary workflows

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

These portable programs produce manuscript PCA, sequence/signature features, statistical comparisons, sampler checks and supplementary panels. They use the [pinned Python environment](../../CONDA_SETUP.md) and audited upstream exports. Current analyses do not require the historical R environment.

## Main analysis and sampler check

After the human, bat COSMIC and experimental-profile tables have been exported and audited, submit from the SeaWulf project root:

```bash
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/author_hpc/validate_manuscript_results_seawulf.slurm
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/author_hpc/validate_signature_sampler_seawulf.slurm
```

The two jobs are independent. Each requests eight CPUs, 8 GB and four hours as allocation limits. Wrappers activate the configured `SINABRO_ENV`; Python entry points accept explicit paths on other computers.

| Program | Calculation | Main output |
| --- | --- | --- |
| `validate_manuscript_results.py` | PCA, GC3, signature features, statistical and species comparisons | `analysis_summary.json`, tables and PDF/PNG panels |
| `validate_signature_sampler.py` | 10,000 random 2,200-bp sequences per each of 79 signatures | `sampler_validation_summary.json`, channel counts, cosine summaries and plots |

The downstream job reads the human export, both bat exports and all three experimental-profile exports. It verifies file manifests, expected populations, signature identity and input compatibility. Default wrapper paths identify the reference campaign; override `SINABRO_PAPER_HUMAN`, `SINABRO_PAPER_BATS`, `SINABRO_PAPER_APOBEC` and `SINABRO_PAPER_OUTPUT` for another project. Output locations must be new.

PCA is mean-centred full SVD without standardization, using the historical sign convention. All 79-, 60- and 41-signature groups are retained; Figure 4 fits the 60-signature model and displays the 41 known-aetiology signatures. The original comparisons include Kruskal–Wallis, median-centred Levene, Dunn/Bonferroni and asymptotic Mann–Whitney tests with the documented family corrections. Transcript-level observations are not independent biological replicates. See [methods and reference results](../../docs/validation_reference.md).

Sampler validation uses 790,000 single-substitution events and an independent conditional probability oracle based on available trinucleotide opportunities. Cosine to the original signature is descriptive; opportunity-adjusted checks verify sampling behavior rather than enforcing the old rounded cosine.

## Supplementary analyses

| Analysis | Execution guide | Interpretation |
| --- | --- | --- |
| S11 paired historical estimates | [S11_HISTORICAL.md](../S11_HISTORICAL.md) | 3,740 matched transcript pairs; archived estimates, not a new full-cohort BLOSUM campaign |
| S12 BRCA1-203 simulation | [BRCA1_GRAPH.md](../BRCA1_GRAPH.md) | 1,000 SBS2/BLOSUM80 paths with reconstructable DNA states and chronological edges |
| COSMIC breast-SNV preparation | [Coordinate preparation](../../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md) | Verified archive/phenotype join and source-to-target genomic coordinate mapping |
| S12 reported-allele annotation | [BRCA1_OBSERVED_SNVS.md](../BRCA1_OBSERVED_SNVS.md) | Single-substitution allele matches; multi-substitution components do not prove observed genotypes |

The reference full campaigns and supplementary stages have passing reports. The [validation reference](../../docs/validation_reference.md) records the numerical results and scope. Run instructions are provided for independent reproduction; already completed compatible campaigns can proceed directly from their audited outputs.

Historical source recovery and old/new comparisons are in [historical/validation](../../historical/validation). Earlier bounded sensitivity examples remain separate from the manuscript's archived S11 population; `run_supplementary_cases.py --help` describes the portable implementation, and the older generic Slurm example is retained under `historical/examples/`.
