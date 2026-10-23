# Reported results

These reports record the authors' completed HPC run (job 2350953). They were
copied from the supplied console output. Script and shared-reader SHA256 values
match the repository; both reports refer to the same previously reported audited
human export. See [provenance.json](provenance.json) for the verification boundary.
Full output CSVs and figures remain on the authors' HPC; their checksums are
listed in the reports. This directory contains report summaries, not the full
gene-by-signature matrices or all per-case precision rows.

## One-step baselines

[one_step_summary.json](one_step_summary.json) reports exact initial synonymous
probabilities for 19,309 fixed gene representatives and 79 signatures, compared
with the existing 1,000-path RUMS estimates. No new trajectories were generated.

| Signature | Uniform-neighbour probability vs RUMS, Spearman rho | Signature-weighted initial probability vs RUMS, Spearman rho |
| --- | ---: | ---: |
| SBS2 | 0.52413 | 0.99509 |
| SBS13 | 0.64033 | 0.98957 |
| SBS19 | 0.47337 | 0.98980 |

Across the 79 signatures, the correlation between the two cohort means
(mean weighted initial probability and mean RUMS) is 0.99886. This is a different
comparison from the within-signature, across-gene correlations in the table.

RUMS closely tracks initial signature-weighted synonymous opportunities. Its
sequential stopping-time definition differs from a first-step probability;
these results do not establish independent biological information or superiority
over the simpler baseline. They support a descriptive, sequence- and
mutation-model-conditioned interpretation. Neither a synonymous-recoding
adaptation null nor fixed-n robustness was calculated.

## Saved-path precision

[stability_summary.json](stability_summary.json) records a metadata-selected
sample of 183 CDSs from the 19,309 representatives, with all 79 signatures
(14,457 combinations). Sampling used length/GC3 tertiles and metadata extrema,
without using RUMS. Existing complete checkpoint paths were checked against
the audit and export before analysis.

| Saved paths used | Median Monte Carlo SE | 90th percentile SE | Largest SE | Median absolute difference from the saved 1,000-path mean | Correlation of sample profile means with 1,000-path profile means |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.06497 | 0.09688 | 0.28628 | 0.041 | 0.99810 |
| 200 | 0.04627 | 0.06784 | 0.19578 | 0.027 | 0.99912 |
| 500 | 0.02954 | 0.04257 | 0.11980 | 0.014 | 0.99954 |
| 1,000 | 0.02087 | 0.03000 | 0.08411 | 0 (by identity) | 1 (by identity) |

SBS2 has the highest sample mean at every prefix; SBS19 is second at every
prefix. At 1,000 paths, SBS2's median SE is 0.03469 and its 90th percentile SE
is 0.05930. These values characterize Monte Carlo variation within this sample,
not variation among independent animals/species or a genome-wide error bound.

The prefixes overlap, and the full saved mean is not independent ground truth.
Agreement with itself at 1,000 paths is an identity. This check does not prove
universal convergence or that 1,000 paths suffice for every sequence/profile.

The two selected computational checks are complete. No PCA scaling comparison
or bat species confidence-interval analysis was performed in this stage.
