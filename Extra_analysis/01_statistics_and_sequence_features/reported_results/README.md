# Verified human signature summaries and manuscript placement

These are the two 79-row descriptive summary tables from the passing full
additional analysis. Their bytes match the SHA256 inventory in the supplied HPC
execution report. `results_summary.json` records the paired-analysis design,
selected results, input/source fingerprints and verified rankings. The full
3,081-row paired table and other outputs remain in the full execution output;
they are not deposited in this small summary directory.

No original manuscript or original-paper analysis file was edited to add this
record. The locations and wording below are proposals for author incorporation.

## Mean and median are distinct summaries

| Cohort | Signature | Mean RUMS | Median RUMS |
|---|---|---:|---:|
| All 73,214 retained transcripts | SBS2 | 0.789106 | 0.679 |
| All 73,214 retained transcripts | SBS19 | 0.719614 | 0.683 |
| One representative CDS per gene, n = 19,309 | SBS2 | 0.833675 | 0.727 |
| One representative CDS per gene, n = 19,309 | SBS19 | 0.755544 | 0.730 |

SBS2 has the highest mean among the 79 profiles in both cohorts. SBS19 has the
highest median in both cohorts. The original highest-mean result is retained;
it should not be rewritten as a highest-median or an unqualified highest-under-
every-summary result. This descriptive ranking alone does not establish
significant paired superiority over every other signature or adaptation.

## Where these results belong

| Location | Material to incorporate |
|---|---|
| Results, paragraph introducing the human signature ranking | Retain the highest-mean SBS2 statement; add a short qualification that SBS19 has the highest median. Distinguish the full-transcript descriptive cohort from gene representatives. |
| Methods, Statistics | Explain longest-retained-CDS selection, matched gene rows, Friedman, 3,081 paired Wilcoxon contrasts, Holm family, paired effects, pointwise median CIs and underflow handling. |
| Supplementary Table S1 replacement package | Include a descriptive section with all 79 profiles and a paired-testing section using `human_gene_paired_comparisons.csv` from the full output. Retain the original all-transcript summary alongside the gene-representative summary. Final subtable numbering is not assigned here. |
| Reviewer 2, point 2.4 | Report the completed repeated-measures correction and quantify effects/intervals. Explain that mean and median rankings differ. |
| Reviewer 2, point 2.2 | Use the direct paired SBS2/SBS13 contrast without generalizing to all APOBEC-associated profiles. |
| Abstract/Highlights | Specify “highest mean RUMS”; the detailed mean-versus-median comparison can remain in Results and the supplementary table. |

## Proposed Results wording

Across the 79 COSMIC SBS profiles, SBS2 had the highest mean RUMS in both the
full transcript cohort (n = 73,214; mean 0.7891) and a cohort retaining one
representative CDS per gene (n = 19,309; mean 0.8337). The median ranking
differed: SBS19 had the highest median in the full cohort (0.683, versus 0.679
for SBS2) and representative cohort (0.730, versus 0.727 for SBS2). Thus, the
highest-mean result is distinguished from the highest-median result.

## Paired-analysis results available for the response

Friedman chi-square was 1,018,044.2265 (78 df; Kendall W = 0.6759464). The
numeric p-value underflowed; it is not an exact zero probability. All 3,081
unordered signature pairs were evaluated with paired Wilcoxon tests and Holm
correction.

For SBS2 versus SBS13, the median of within-gene differences was 0.503
(pointwise 95% median interval 0.498–0.509); matched rank-biserial effect was
0.99999976. Marginal medians were 0.727 and 0.212, respectively. These results
concern the operational protein-sequence criterion and conditional sequence/gene
sampling; they do not measure functional loss or prove evolved protection.

The report does not supply a verified gene-based CI for bat species contrasts.
Those remain a separate unresolved part of the requested supplementary-table
expansion. The full pairing table is required before drawing conclusions about
SBS2 versus SBS19 beyond the descriptive mean/median values above.
