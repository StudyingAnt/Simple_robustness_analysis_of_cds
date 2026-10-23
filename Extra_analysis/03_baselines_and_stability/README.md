# One-step baselines and saved-path stability

These two checks complete the selected additional computational scope. Both read
existing inputs/results and generate no new RUMS trajectories. The existing
package, original-paper workflows and their outputs are not edited.

Use the repository's pinned Conda environment. From the repository root:

```bash
python Extra_analysis/03_baselines_and_stability/analyze_baselines.py \
  --human-export /PATH/TO/audited_human_export \
  --signatures /PATH/TO/COSMIC_v3.3.1_SBS_GRCh37.txt \
  --out-dir /PATH/TO/new_one_step_results

python Extra_analysis/03_baselines_and_stability/analyze_stability.py \
  --human-export /PATH/TO/audited_human_export \
  --out-dir /PATH/TO/new_stability_results
```

The export directory must retain every file in its authenticated output
inventory. `--human-fasta` can point to a relocated, hash-identical CDS-plus FASTA.
The stability check additionally requires the original saved run and full audit;
their paths are read from the export manifest, or supplied with `--run-dir` and
`--audit-dir`. Checksums and configuration must still match. Outputs must be new
directories outside input directories.

Default production checks require 73,214 human transcripts, 19,309 fixed gene
representatives, 79 signatures and 1,000 trajectories. Explicit software validation
can use `--allow-subset --n-traj N`; stability also needs valid `--prefixes` ending
at N. These bounded runs are not full-paper results.

## 1. Exact first-step baselines

The analysis uses the same longest-retained-CDS gene representatives as the
paired-statistics analysis, with accession tie breaking. Selection does not use
RUMS values.

For each CDS of length L nucleotides, all 3L single-base substitutions are
classified using the standard genetic code and the protein-identity criterion.
The uniform baseline is the fraction that preserve the encoded amino-acid
sequence. Terminal and internal stop symbols are retained, consistent with the
saved simulation criterion; the two flanking bases provide context but are not
mutation targets.

The signature-weighted baseline is the total initial synonymous event weight
divided by all available initial event weight. SBS96 channels are expanded with
equal reverse-complement factors, which cancel in the ratio. Opportunity counts
are computed exactly, without Monte Carlo mutation sampling.

Outputs include:

- `uniform_one_step_baseline.csv`: counts and uniform synonymous fraction.
- `weighted_one_step_baseline.csv`: representative CDS by signature probability.
- `baseline_correlations_by_signature.csv`: Spearman associations with saved
  RUMS, absolute rank shifts and tie-inclusive top-5% overlap counts.
- `rank_discordant_examples.csv`: ten largest rank discrepancies per profile.
- `one_step_baseline_comparison.pdf/.png` and `analysis_summary.json`.

This is a quantitative one-step sequence baseline inspired by neutral-neighbour
robustness definitions. A first-step probability and a sequential stopping-time
expectation are different quantities. The calculation is not fixed-n robustness,
a protein-function assay, an adaptation null, a synonymous-recoding experiment,
or multivariable adjustment for codon usage/context. Correlations are descriptive;
they do not establish an independent causal contribution of RUMS.

## 2. Saved-path precision and prefix sensitivity

The sample is selected from the 19,309 representatives using length and GC3
tertiles. Up to 20 representatives are chosen from each of nine bins by a
deterministic SHA256 ordering with seed 42. Metadata extremes (minimum/maximum
length and GC3) are also included, yielding at most 184 sequences. All 79 stored
profiles are used. RUMS values never enter sample selection.

For each selected checkpoint, its bytes must match the full audit inventory and
its complete path counts must reproduce the exported means. The prefixes of
100, 200, 500 and 1,000 saved paths provide:

- Mean RUMS, path SD and empirical Monte Carlo SE = SD / sqrt(n).
- Proportion of individual paths with zero accepted synonymous steps. This is
  distinct from the fraction of transcript means that equal zero.
- Differences from the full saved 1,000-path mean.
- Profile-mean rank agreement within the selected sequence sample.

Outputs include `path_precision_by_case_and_prefix.csv`,
`precision_summary_by_prefix.csv`, `precision_by_signature_and_prefix.csv`,
`sampled_representatives.csv`, a precision figure and `analysis_summary.json`.

Nested prefixes share observations. The 1,000-path mean is not an independent
truth estimate, and the zero difference at the full prefix is an identity, not
evidence of convergence. Empirical SE may fluctuate between prefixes; no
monotonic decrease is imposed. Zero sample variance does not prove zero true
variance. This is representative numerical precision reporting, not a claim
that 1,000 paths universally suffice for every gene/profile. Gene/species sampling
uncertainty is distinct from Monte Carlo variation.

## Verification and author execution

```bash
python -m unittest discover -s Extra_analysis/03_baselines_and_stability/tests -v
```

Tests enumerate every single-base neighbour of all 64 codons, compare weighted
probabilities with the independent package event matrix, check prefix formulas,
and exercise both workflows using synthetic authenticated exports. Deterministic
metadata-only selection and tampered checkpoint/audit rejection are checked.

See [author_hpc/README.md](author_hpc/README.md) for the authors' SeaWulf submission
example. Other users may run the portable commands directly or adapt scheduler
settings to their own system. No PCA scaling analysis or bat confidence-interval
analysis is included in this stage.
