# Original manuscript validation

These programs reproduce downstream analyses already present in the human/bat manuscript. They use the pinned Python/NumPy/Biopython/SciPy environment. They do not require R, pandas, scikit-learn or a new environment installation. The package and all completed whole-genome trajectories are unchanged.

## Submit two independent jobs

From the SeaWulf project root, after updating this analysis branch:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/validate_manuscript_results_seawulf.slurm
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/validate_signature_sampler_seawulf.slurm
```

The two jobs can run concurrently. Each uses eight shared CPUs; requested time limits are four hours, not measured runtimes. They activate the existing PROJECT/sinabro-revision prefix. This directory name is only an environment location; it does not change the analysis definitions. No previously completed trajectories are rerun.

## PCA, signature features, original statistics and species comparisons

`validate_manuscript_results.py` reads the audited human export (`runs/aggregation_2305071`), bat COSMIC exports (`runs/bat_cosmic_export_2328086`) and experimental-enzyme export (`runs/experimental_apobec_export_2338741`). Every exported input file is checked against its manifest. Full cohort counts, signature-file identity, DNA identity across comparisons and original profile names/order are checked. `--allow-subset` is only for explicitly labelled bounded software checks.

Output defaults to `results/manuscript_validation_JOBID/`. It contains:

- PCA variance/loadings for all 79, no-artifact 60 and known-aetiology 41 signature groups; PC1–PC10 transcript scores, GC3 and whole-CDS GC features.
- Signature entropy (with zero-weight channels handled explicitly), GC targeting and transition/transversion fractions; original transition-versus-mean regression and entropy/SD correlations.
- Kruskal–Wallis, median-centered Levene and two-sided Dunn/Bonferroni comparisons for each signature group.
- Pairwise asymptotic Mann–Whitney/Wilcoxon rank-sum comparisons, with tie correction and continuity correction, for shared COSMIC controls and all six experimental enzyme profiles. Both raw and within-comparison-family Bonferroni p-values are saved; original table adjustment conventions are not assumed to have been recovered.
- The native-process Figure 7g comparison uses **human COSMIC SBS2**, Molossus's own APOBEC3A and Myotis's own APOBEC3A3. It is distinct from applying one fixed experimental enzyme to all three genomes. Mean, SD, median, zero fraction, two-sided results and the original one-sided comparison directions are exported without enforcing the old conclusion.
- PDF/PNG panels for Figures 3–5, bat COSMIC controls, the native-process panel and six experimental enzyme comparisons. All observations and distribution tails are retained. These are regenerated scientific plots, not pixel-identical copies of historical ggplot2 panels.

PCA follows the original default scikit-learn design: mean-centered full SVD with no standardization, using the same svd_flip sign convention. The implementation uses NumPy so the checked environment needs no extra dependency. First-ten scores are exported because those are the components analyzed in the paper. All loadings and variance ratios are retained.

These tests reconstruct the original unpaired transcript-level analyses. Repeated measures across signatures and transcript isoforms do not constitute independent biological replicates; the output explicitly records this limitation. They are not a substitute for a revised repeated-measures or orthologue-matched design. A floating-point p-value of zero is numerical underflow, not an exactly zero probability. Dunn implements pooled ranks with tie correction and two-sided normal tails, followed by Bonferroni; see the [R dunn.test implementation](https://github.com/cran/dunn.test/blob/master/R/dunn.test.R). Median-centered Levene uses the [SciPy 1.8.1 definition](https://docs.scipy.org/doc/scipy-1.8.1/reference/generated/scipy.stats.levene.html).

Review `analysis_summary.json`: `analysis_passed: true` means input checks and computation completed. Actual PC variance, regression coefficients and native-panel results must be read before claiming agreement with the manuscript. Historical assertions are not used as pass conditions. The program writes a new directory atomically and refuses existing output locations.

Override input locations with `SINABRO_PAPER_HUMAN`, `SINABRO_PAPER_BATS` and `SINABRO_PAPER_APOBEC`, and output with `SINABRO_PAPER_OUTPUT`. On other computers run `validate_manuscript_results.py --help` and provide those paths directly.

## Original Figures S1/S2 sampler validation

`validate_signature_sampler.py` generates 10,000 independent 2,200-bp random sequences with expected GC fraction 0.4 for each of the 79 original signatures, and calls the actual Sinabro signature sampler once on each sequence. Total: 790,000 single substitutions. The source package/profile file and seeds are fingerprinted. Jobs use spawned workers with independent signature-keyed streams.

Output defaults to `results/signature_sampler_JOBID/`. It includes channel counts/expected counts, per-signature cosine similarities, 79 SBS96 profile overlays and the cosine overview. This is a separate original-paper validation experiment, not another genome-wide RUMS calculation.

The observed fractions are compared with **an independent conditional oracle**, computed by enumerating forward/reverse-complement trinucleotide opportunities and their original profile weights. Original-profile cosine similarity is descriptive because realized fractions depend on available sequence contexts. The historical average cosine of 0.974 is not enforced. Common-bin residual checks and conservative Poisson-binomial Chernoff bounds check both common and rare channels; the threshold accounts for 79 × 96 channels. Source-free impossible events fail explicitly. `sampler_validation_summary.json` must show `validation_passed: true` before citing sampler validation as passed.

Both programs use raw research inputs and portable Python code; cluster settings are isolated in Slurm files. The scientific package is imported by the sampler job only. The downstream job reads data and performs no Sinabro simulation.

## Remaining original supplementary analyses

S11 (BLOSUM80/SBS2 subset) and S12 (BRCA1 genotype/phenotype graph) are original-paper analyses. They remain separate from the two jobs above. **Archive recovery update:** the [historical supplementary input audit](HISTORICAL_SUPPLEMENTARY_INPUTS.md) identifies 3,740 S11 results with their exact transcript sequences, and S12 graph code/artwork. The COSMIC observed-variant overlay table and matching rule remain unresolved. Prioritize these recovered sources for original-panel validation; the examples below remain separate bounded reconstructions. A regenerated analysis must state its chosen subset/settings and data provenance rather than claiming byte-identical historical reproduction. S12 red observed-variant points cannot be reconstructed without the source variant data. These panels must not be marked verified solely because the package's BLOSUM unit tests passed.

The following optional examples were prepared before the archive was recovered. They can run after the two core validations, but do not replace validation of the recovered historical sources:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
SINABRO_SUPPLEMENT=s11 sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/run_supplementary_cases_seawulf.slurm
SINABRO_SUPPLEMENT=s12 sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/run_supplementary_cases_seawulf.slurm
```

Each job requests eight CPUs, 8 GB and four hours. These are allocation limits, not measured runtime estimates. The wrapper uses existing inputs and the existing Conda environment. It does not rerun whole-genome trajectories. Portable command-line usage is available through `run_supplementary_cases.py --help`.

Both reconstructions use SBS2 and a BLOSUM80 global protein alignment, gap-open −10, gap-extension −0.5 and stopping distance ≥20. The path cap is 10,000 proposed mutations; reaching it without a stopping event fails the run rather than contributing a truncated robustness value. The `score` stop policy evaluates translated stop-containing sequences through the alignment score. BLOSUM distance is a proxy and is not an experimental measurement of functional loss.

- **S11:** the historical exploratory selector yields 985 CGC-v99 gene-201/gene-203 transcripts. This is a different exploratory source from the subsequently recovered 3,740-result Figure S11 input. The frozen [case manifest](supplementary_cases.json) chooses a new, explicit set of 12 length-quantile representatives (291–7,029 coding nucleotides). Each receives 1,000 paired-seed paths under protein identity and under the BLOSUM criterion: 24,000 paths total. Shared mutation prefixes and the relaxed-path length bound are checked. The output includes per-path counts, criterion means, paired-difference Monte Carlo standard errors and PDF/PNG scatter plots. This small subset illustrates criterion sensitivity; it does not establish a genome-wide effect.
- **S12:** BRCA1-203 (ENST00000357654.9; 5,592 coding nucleotides) receives 100 BLOSUM paths. Graph nodes represent actual simulated DNA sequences, with DNA substitutions, Hamming distances and BLOSUM distances exported alongside edges. The display uses a deterministic radial layout, not the historical exploratory t-SNE layout. The graph contains **simulation only**; no observed cancer variants are imputed or substituted from another database.

Read `results/supplementary_s11_JOBID/supplementary_summary.json` and the corresponding `s12` file. `simulation_passed: true` confirms completed checked paths. `original_figure_reproduced: false` records that these new examples are not the original figure reconstruction; it remains false even when an example succeeds. For S12, `observed_cosmic_variant_overlay_reproduced` also remains false. New output directories are required and input hashes are checked before and after computation.
