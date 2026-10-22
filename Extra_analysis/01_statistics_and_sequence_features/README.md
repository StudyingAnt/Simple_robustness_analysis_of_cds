# Paired statistics and sequence features

This analysis reuses the audited RUMS exports. It generates no new trajectories
and changes no original-paper code, input or output. Use the pinned environment
in the repository's `CONDA_SETUP.md`; no additional dependencies are required.

## Run on any computer

From the analysis repository root, with the environment active:

```bash
python Extra_analysis/01_statistics_and_sequence_features/run_analysis.py \
  --human-export /PATH/TO/human_export \
  --bat-export /PATH/TO/bat_cosmic_export \
  --apobec-export /PATH/TO/experimental_apobec_export \
  --out-dir /PATH/TO/new_additional_analysis_output
```

`human_export` contains `aggregation_manifest.json`,
`rums_mean_by_transcript_id.csv`, `transcript_annotations.csv` and the other
authenticated exporter outputs. `bat_cosmic_export` contains `molMol2/` and
`myoMyo6/` exports; `experimental_apobec_export` contains all three species.
Keep all files listed in each manifest, since their hashes are verified.
The defaults require 1,000 paths per combination, 79 human COSMIC signatures,
eight bat COSMIC signatures, six enzyme profiles and the full frozen populations.

The manifest-recorded FASTAs must be present. If files were moved, supply
`--human-fasta /PATH/TO/human.cdsplus.fa`, `--molmol2-fasta /PATH/TO/molMol2.cdsplus.fa`
and `--myomyo6-fasta /PATH/TO/myoMyo6.cdsplus.fa`; their hashes must match the
simulation inputs. Changing a path does not permit changing a sequence.

For a human-only first run, omit both `--bat-export` and `--apobec-export`.
For an explicitly bounded software check, add `--allow-subset --n-traj N` with
the actual path count. Such a run is marked as a subset and is not full-paper
inference. Existing output directories are refused.

## Statistical design fixed before inspecting the additional results

- Human inference uses one representative transcript per stable gene ID: the
  longest retained CDS, with versioned transcript accession as a lexical tie-break.
  Selection never uses RUMS. The same representative rows are used for every
  signature. Full-transcript summaries remain available separately.
- Friedman is the tie-corrected omnibus test across all 79 signatures. Paired
  two-sided Wilcoxon signed-rank contrasts cover all 3,081 unordered pairs;
  Holm correction uses this entire family, not a selected subset of favourable
  comparisons. The analysis uses the normal approximation, drops zero
  differences, averages tied ranks and applies no continuity correction.
- The signed-rank null concerns the distribution/symmetry of differences; it is
  not a general test solely of equal medians. Effect magnitude and paired median
  differences are therefore reported explicitly rather than inferred from p alone.
- Means from complete paths lie on the `1 / n_trajectories` lattice. Integer
  path-count sums recover this lattice before ranking differences, preventing
  floating-point subtraction from breaking ties. Off-lattice inputs are refused.
- Positive matched rank-biserial effect and positive paired median difference
  both mean `group1 > group2`. The paired median is `median(x - y)`, which is
  generally different from `median(x) - median(y)`.
- Median CIs use binomial/order-statistic limits, retaining zero differences.
  They are pointwise 95% intervals, not simultaneous intervals for all pairs and
  not Wilcoxon/Hodges–Lehmann location intervals. Discreteness can make coverage
  conservative; the achieved continuous-case coverage bound is exported. For
  very small validation subsets, an unbounded limit is recorded as a blank/null.
- These are sequence/gene sampling summaries, conditional on the retained
  annotation and simulated estimates. One CDS per gene reduces isoform
  duplication; it does not create independent animal samples, remove all
  evolutionary dependence, or propagate per-path Monte Carlo uncertainty.
- Tiny paired p-values retain log10 p and log-space Holm values. A numeric zero
  is explicitly marked as underflow. The omnibus p is also flagged on underflow;
  no exact zero probability is claimed.

## Composition, species comparisons and outliers

GC3 is calculated from the third bases of the CDS after excluding the single
flanking nucleotide at each end. The terminal stop codon is retained, consistently
with the frozen CDS inputs; internal-stop transcripts remain present and flagged.
Length, total GC and GC3 summaries/ECDFs describe the retained transcript cohorts;
they are not an orthologue matching or covariate-adjustment analysis.

For each human profile, the default outlier file lists the ten largest RUMS
means, with transcript, gene, length, GC3 and internal-stop flag. Ties are ordered
by accession. `--top-outliers` changes this display count, not the inferential
cohort. Spearman associations with GC3 and length are descriptive and have no
added significance tests.

Species tables distinguish all eight shared fixed COSMIC profiles, six fixed
experimental profiles, and the native diagonal that changes both cohort and
profile. Overlapping table entries can be taken from this single canonical file.
They contain medians/IQRs, differences of cohort medians and independent
rank-biserial effects. No paired test is applied to unmatched species.

The current frozen bat annotations do not document gene grouping. The species
tables therefore have no inferential p-values or CIs; these are explicitly blank,
not zero. This does **not** fulfil the request for species-table CIs. Verified
grouping would be needed for gene-based sampling intervals, whose interpretation
would still be distinct from animal/species replication. No claim of coevolution,
native-process adaptation or lifespan causation follows from these comparisons.

## Outputs and reviewer coverage

| File | Purpose |
|---|---|
| `additional_analysis_summary.json` | Completion status, Friedman result, input/source hashes, cohort sizes, design and limitations |
| `human_representative_transcripts.csv` | Reproducible gene-to-CDS selection |
| `human_gene_paired_comparisons.csv` | S1 replacement candidate: paired tests, Holm p, effects and median CIs |
| `human_gene_signature_summary.csv` | Per-signature gene-cohort medians and distributions |
| `human_transcript_signature_summary.csv` | Original full-transcript descriptive summaries |
| `*_sequence_features.csv`, `sequence_feature_summary.csv` | Three-cohort composition and length |
| `gc3_cohort_distributions.pdf/.png` | Full-cohort GC3 empirical distributions |
| `human_top_outliers_by_signature.csv` | Figure S3 annotation candidates |
| `human_composition_correlations.csv` | Descriptive GC3/length associations |
| `species_profile_distributions.csv`, `species_descriptive_contrasts.csv` | Canonical descriptive S3/S4 candidates; only in the three-species run |

These outputs address Reviewer 2's repeated-measures concern (2.4), composition
limitation (3), outlier request (6), and existing SBS2/SBS13 contrast (2.2).
Effects, uncertainty and composition summaries can also support overlapping
Reviewer 1 concerns. They do not supply a synonymous-recoding null, a numerical
comparison with earlier robustness measures/simulators, a species-by-process
orthologue model, assay replicate metadata, or a phylogenetic lifespan test.

## Verification and author execution

```bash
python -m unittest discover \
  -s Extra_analysis/01_statistics_and_sequence_features/tests -v
```

Tests compare signed ranks to SciPy, check exact binomial coverage, verify effect
orientation and Holm adjustment/underflow, and exercise the complete workflow
with authenticated synthetic exports. Synthetic validation is distinct from the
authors' production results. The authors ran the main calculations on SeaWulf;
see [the author-only submission example](author_hpc/README.md).

Method documentation:
[SciPy Friedman](https://docs.scipy.org/doc/scipy-1.8.1/reference/generated/scipy.stats.friedmanchisquare.html),
[SciPy Wilcoxon](https://docs.scipy.org/doc/scipy-1.8.1/reference/generated/scipy.stats.wilcoxon.html),
[R Holm adjustment](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html),
[Geyer: binomial/order-statistic median intervals](https://www.stat.umn.edu/geyer/old/5601/examp/sign.html).
