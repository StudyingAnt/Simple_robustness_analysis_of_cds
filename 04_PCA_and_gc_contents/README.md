# PCA, coding-sequence features and signature features

The current portable implementation is [validate_manuscript_results.py](../07_manuscript_validation/validate_manuscript_results.py). Follow the [manuscript analysis guide](../07_manuscript_validation/README.md) for command-line and SeaWulf usage. It reuses audited exported tables and performs no new whole-genome simulation.

The analysis retains the original 79-signature group, the 60-signature group without sequencing artifacts, and the 41-signature group with known aetiology. Group identities and ordering are frozen in [manuscript_signature_groups.json](../02_run_sinabro/manuscript_signature_groups.json).

PCA uses mean-centered full SVD without feature standardization, with the original sign convention. All loadings and explained-variance ratios are exported, together with the first ten transcript scores, GC3 and whole-CDS GC features. Entropy, transition fractions and the corresponding regressions are computed from the pinned signature table. Actual observed coefficients and variance ratios are recorded rather than enforced to match the manuscript.

The Python and R files in this directory preserve the historical analysis. They may depend on fixed paths, historical tables and a separate R environment. The combined portable implementation provides the current reproduction entry point; the old files remain available for provenance and method comparison.
