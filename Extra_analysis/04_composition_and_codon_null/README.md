# Sequence composition and a synonymous-codon null

These are additional human-CDS analyses. They keep the original paper workflows,
saved outputs and Sinabro implementation unchanged. One analysis jointly models
sequence covariates using existing RUMS estimates; the other generates new
trajectories after protein-preserving synonymous-codon resampling. They answer
different questions.

Completed numerical results, figures and all additional null path counts are in
[`reported_results/README.md`](reported_results/README.md). The reported null
effects and confidence intervals can be checked directly from that bundle using
`verify_reported_results.py`, without rerunning the simulations.

Use the repository's pinned environment and an installed Sinabro 3.0.0 checkout.
Python can run directly on an ordinary multicore computer. Eight worker processes
were used by the authors; no institutional scheduler is required.

The portable commands below are the entry points for reuse. The authors' local
command sequence is also kept separately in
[`author_execution/run_local.sh`](author_execution/run_local.sh). Activate the
environment and supply the four input/output path variables shown in that file.

## Prepare authenticated inputs

From the repository root:

```bash
python Extra_analysis/04_composition_and_codon_null/prepare_inputs.py \
  --human-export /PATH/TO/audited_human_export \
  --human-fasta /PATH/TO/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --signatures /PATH/TO/sinabro/sinabro/data/COSMIC_v3.3.1_SBS_GRCh37.txt \
  --out-dir /PATH/TO/new_composition_inputs
```

Every exported file is checked against the aggregation manifest. The CDS-plus
FASTA and signature table must match the original simulation checksums. The
original full export is required: 73,214 transcripts, 79 profiles and 1,000 paths
per combination. The fixed representative for each of 19,309 genes is its longest
retained CDS, breaking length ties by transcript accession. RUMS is not used in
selection.

Five profiles are fixed before new simulations: SBS2 (focal APOBEC C>T profile),
SBS13 (APOBEC-associated C>G comparison), SBS19 (highest median in the existing
representative-gene results), SBS1 (another transition-dominated profile) and
SBS4 (transversion-rich comparison). This is not an all-profile null analysis.

## 1. Joint composition analysis

```bash
python Extra_analysis/04_composition_and_codon_null/analyze_composition.py \
  --inputs /PATH/TO/new_composition_inputs \
  --out-dir /PATH/TO/new_composition_results
```

The responses are the existing mean RUMS for each gene/profile. The primary scale
is `log1p(RUMS)`; raw RUMS is a sensitivity analysis. Nested feature sets are:

1. Log CDS nucleotide length, total GC and GC3.
2. Those features plus the relative frequencies of all 64 codons.
3. Those features plus the relative frequencies of all 64 trinucleotide contexts
   centered on mutable CDS bases, including context supplied by the flanks.
4. Those features plus the exact initial signature-weighted synonymous probability.

Frequency features represent codon usage/composition and initial context
opportunities, not expression-adjusted translation efficiency or experimental
mutation rates. The codon frequencies also contain amino-acid composition.

GC, GC3 and frequency features are structurally dependent. Standardized ordinary
least squares therefore uses SVD rank reduction with a relative singular-value
cutoff of `1e-10`. No individual coefficient, coefficient significance or
independent causal contribution is interpreted. Effective ranks are reported.

Five fixed SHA256 gene folds are shared across profiles. Scaling and regression
are fitted on training genes only. The main fit metric is pooled gene-held-out
R-squared, distinct from training R-squared. A separate partial-rank analysis
correlates residual ranks of RUMS and initial synonymous probability after both
are projected on the ranked composition features. It is a conditional
association, not a causal decomposition or evidence of adaptation.

Outputs include model metrics, held-out predictions, partial-rank associations,
`composition_adjustment.pdf/.png`, and `analysis_summary.json`. Analyses are
within the human cohort; they do not standardize species or orthologues.

## 2. Protein-preserving synonymous-codon null

```bash
python Extra_analysis/04_composition_and_codon_null/run_codon_null.py \
  --inputs /PATH/TO/new_composition_inputs \
  --sinabro-source /PATH/TO/sinabro \
  --out-dir /PATH/TO/new_codon_null_simulation \
  --recodings 20 --n-traj 500 --workers 8

python Extra_analysis/04_composition_and_codon_null/analyze_codon_null.py \
  --inputs /PATH/TO/new_composition_inputs \
  --simulation /PATH/TO/new_codon_null_simulation \
  --out-dir /PATH/TO/new_codon_null_results
```

Six representatives are sampled without replacement in each of 16 length-quartile
by GC3-quartile strata, giving 96 genes. The selection seed is 20261025 and does
not use RUMS. The null samples each codon uniformly and independently from the
synonymous codons for its amino-acid symbol, using the standard genetic code.
Stop symbols are retained and their synonymous stop codons can be resampled.
Protein sequence, CDS length and both flanking nucleotides are preserved.
GC, GC3, codon frequencies and context frequencies are allowed to change.

Twenty recoded sequences and the original sequence each receive 500 independent
paths for each of five profiles: 5,040,000 new paths. This additional experiment
uses 500 paths per sequence/profile; the existing original-paper calculation
continues to use 1,000. Each path uses the unchanged Sinabro protein-identity
criterion and updates mutation context after substitutions. The bounds are
`cds_start=1`, `cds_end=len(sequence)-1`; `max_steps=10000` is a safety cap, not
an imputed trajectory length. Any incomplete path fails the run.

Recoding and path RNG streams are independently indexed by master seed,
transcript, recoding number and profile. Recoding variants use the same starting
DNA across profiles. The complete trajectories' accepted-step counts, exact
first-step probabilities, sequence/protein hashes and configuration are saved.
The source package must match the authenticated original calculation.

The effect is original mean RUMS minus the mean over recoded sequences. Genes,
not recodings or paths, are the units for the 10,000 stratified gene-bootstrap
replicates. Stratum estimates are weighted by their proportions in the full
19,309-gene cohort. This yields a limited sampled-cohort estimate; it does not
replace a full-cohort null. Confidence intervals treat finite simulation
estimates as the gene-level observations. A separate conditional Monte Carlo SE
reports original-path and finite-recoding/path uncertainty.

The SBS2 excess is paired against the excess for each comparison profile.
Pointwise 95% intervals and Bonferroni-adjusted intervals across the four
comparisons are reported. A first-step probability audit checks empirical
`P(accepted steps >= 1)` against exact sequence-specific probabilities. If a
timing run is supplied with `--timing-dir`, shared variants must reproduce their
saved path prefixes exactly.

An excess over this uniform-codon null alone does not identify natural selection
or establish evolved protection. It may reflect codon or context composition
that this null deliberately changes. The gene/profile coverage, sampling
uncertainty and finite simulation budget must accompany interpretation.

## Verification

```bash
SINABRO_SOURCE=/PATH/TO/sinabro python -m unittest discover \
  -s Extra_analysis/04_composition_and_codon_null/tests -v
```

Checks cover a known regression response with redundant covariates, isolation of
held-out outcomes, preserved protein including stop symbols, exact stochastic
replay, nested RNG prefixes, and agreement with exact first-step probabilities.
