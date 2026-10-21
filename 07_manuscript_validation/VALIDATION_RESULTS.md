# Full manuscript downstream and sampler validation

The author supplied both complete HPC JSON summaries on 2026-10-21. Both validations passed. Their script and signature/group-file checksums match the current checked sources. The downstream report is a full-data run (`bounded_subset: false`), not a pilot. This records the numerical evidence; generated HPC artwork and every exported CSV have not been independently opened or rehashed locally from these JSON attachments.

## Signature sampler: Figures S1/S2

- All 79 original signatures, each with 10,000 independent 2,200-bp sequences at expected GC 0.4, one substitution per sequence: 790,000 substitutions total.
- All 79 independent conditional-oracle checks passed; no failed signatures.
- Mean cosine similarity to the original COSMIC profiles: **0.9751868632**, versus **0.974** reported in the manuscript. The historical exact random streams are unavailable.
- Mean cosine similarity to the independently enumerated sequence-conditional expectation: **0.9992583824**; minimum **0.9970532197**.
- The lowest original-profile cosine is SBS94 (**0.9229309498**); the highest is SBS10a (**0.9995766633**). Update historical minimum/maximum values and their labels when replacing the old numerical figure annotations.

The two cosine comparisons answer different questions. Original-profile similarity is descriptive; sequence-context opportunity affects realized mutation fractions. Sampler correctness was checked against the conditional oracle rather than requiring the old mean or a fixed cosine threshold.

## Human PCA and signature-feature analyses

The principal manuscript comparison uses the 60-signature no-artifact group. The historical reference below is an independent recomputation from the preserved historical matrix.

| Quantity | Historical matrix | New full HPC calculation |
|---|---:|---:|
| PC1 + PC2 explained variance | 73.24703% | 73.23996% |
| PC1–PC10 explained variance | 90.20846% | 90.20196% |
| Transition fraction versus mean RUMS, R² | 0.71495226 | 0.71485790 |

Thus the manuscript's rounded **73%**, **90%** and **R² = 0.715** remain supported. The new transition regression has slope 0.34949867 and p = 1.92833 × 10⁻¹⁷. PC1/GC3 correlations remain positive (Pearson 0.89825; Spearman 0.92443). Entropy/SD correlations remain negative (Pearson −0.68647; Spearman −0.75688).

Kruskal–Wallis and median-centered Levene calculations completed for the 79-, 60- and 41-signature groups. Reported floating-point p-values of zero are underflow, not exact zero probabilities. These are original transcript-level analyses, with the independence limitations already documented in the execution guide.

### Figure 4 loading-table comparison completed

The complete new `pca_noseqerr_loadings.csv` (3,600 coefficients) was supplied and its SHA256 exactly matches the output entry in the earlier downstream report: `c958c30b515b562a69b1b87bad09352c529295ee722e7f6d17b687a6416bf7a7`. All 60 signatures occur once per component and the loading matrix passes an orthonormality check. Comparison uses the original Figure 4 rule: fit the 60-signature PCA and display only the 41 known-aetiology signatures.

| PC | Historical displayed positive maximum | New displayed positive maximum |
|---|---|---|
| PC1 | SBS2 | SBS2 |
| PC2 | SBS2 | SBS2 |
| PC3 | SBS1 | SBS1 |
| PC4 | SBS1 | SBS1 |
| PC5 | SBS10a | SBS10b |
| PC6 | SBS7d | SBS7d |
| PC7 | SBS7a | SBS7a |
| PC8 | SBS84 | SBS84 |
| PC9 | SBS22 | SBS22 |
| PC10 | SBS1 | SBS1 |

Only PC5's displayed maximum changes. Its two leading coefficients remain close:

| Signature / PC5 | Historical coefficient | New coefficient |
|---|---:|---:|
| SBS10a | 0.26305612 | 0.26022800 |
| SBS10b | 0.26055271 | 0.26070322 |

The new lead is 0.00047522. For each of PC1–PC10, the full 60-coefficient historical/new vector cosine is at least **0.99923324**, with the same signs under the recorded convention. This supports very similar loading patterns, while lower ranks can still switch; it does not establish which implementation or random-stream change caused each difference.

Update the Figure 4 text/caption to identify SBS10b as the new PC5 maximum and note the similarly large SBS10a loading. PC8/SBS84, PC9/SBS22 and PC10/SBS1 remain the displayed maxima; the earlier JSON's different global maxima include signatures excluded from the Figure 4 display. The numerical loading check is complete; visual inspection of the new final artwork remains separate.

## Native-process species panel

This panel compares human COSMIC SBS2, Molossus's own APOBEC3A and Myotis's own APOBEC3A3. It is not a fixed-profile, orthologue-matched comparison.

| Dataset | Transcripts | Mean RUMS | Median RUMS |
|---|---:|---:|---:|
| Human | 73,214 | 0.789106 | 0.679 |
| Molossus molossus | 31,185 | 0.665150 | 0.648 |
| Myotis myotis | 23,213 | 0.787936 | 0.727 |

The original rank-test directions are retained: human > Molossus, Myotis > human, and Myotis > Molossus. Human/Myotis two-sided p is 4.69767 × 10⁻⁵⁶ (Bonferroni 1.40930 × 10⁻⁵⁵), with rank-biserial effect −0.06860 in the human-minus-Myotis orientation. This is a small rank effect despite the small p-value and large transcript count. The test is not a test that one arithmetic mean is larger, nor proof of biological adaptation.

Myotis has the highest median, but its mean is **0.00117019 lower than the human mean** in this native-process panel. The text must distinguish mean, median and rank-distribution comparisons; an unqualified claim that Myotis has the highest mean would be incorrect for these results.

## Historical S11 matched-population analysis

Full HPC job **2345768** passed: all **3,740** input sequences and paired transcript identities verified. The source data, manifest and script hashes match the local verified source. Both criterion distributions now use the same 3,740 transcripts, compared with 73,214 versus 3,740 in the historical plotted populations. Historical means are 0.8116874332 (identity) and 3.8720010695 (BLOSUM); medians are 0.705 and 3.594. All recovered BLOSUM estimates are at least their corresponding identity estimates.

This is reproducible postprocessing of archived simulation estimates. `fresh_simulation` and `corrected_package_blosum_recalculation` remain false by design; those flags do not indicate a failed analysis. Historical per-path seeds were not recovered. See [S11 execution and scope](S11_HISTORICAL.md). The author supplied the successful console summary; HPC-generated figures were not independently inspected in this verification.

## Provenance and remaining scope

- Shared signature table SHA256: `9d9c212cc14aa660da1b7fd9e60e696c8e2c885a7606c27d20ad8997954554ea`.
- Human input SHA256: `c2f30678cb9f2252bd4717d1f0b8fef9e9575b4be843aa64478c4b1e0422a751`.
- Sampler source SHA256: `2413daa2187f8012de881c4a790373d700ccb5806844f5eda6a1d9459934fcf6`.
- Downstream source SHA256: `c4847612066dc84ccd21d362e08a7338b07fbf82223ee59447a9d1eabde05533`.
- Source reports list 162 sampler output files and 31 downstream output files. The supplied reports do not include the Slurm job IDs or their own absolute output directories.

The downstream JSON's static `remaining_original_analyses` list includes the separate S1/S2 job. That list is not a cross-job status check: S1/S2 is now verified by its own successful report and must not be rerun merely because it appears there.

The Figure 4 numerical loading review is complete. Finish final artwork checks and address the [recovered S11 sources and unresolved S12 observed-variant mapping](HISTORICAL_SUPPLEMENTARY_INPUTS.md). These successful core reports do not by themselves establish complete supplementary reproduction or a final tagged paper release.
