# Methods and reference validation

This record summarizes the checked manuscript workflows and their scientific scope. Reference-campaign counts and numerical results come from preserved HPC manifests and author-supplied complete reports; large checkpoint stores and generated figures are maintained outside Git. Source/input hashes connect those reports to the checked code. Detailed development comparisons remain in the [historical validation record](../historical/validation/DEVELOPMENT_VALIDATION_RECORD.md).

## Software and inputs

The reference environment is Python 3.9.15, NumPy 1.21.6, Biopython 1.79, SciPy 1.8.1 and Matplotlib 3.6.2 with Sinabro 3.0.0. [software_manifest.json](../software_manifest.json) records the scientific package source fingerprint and code reference. Human CDS-plus SHA256 is `c2f30678cb9f2252bd4717d1f0b8fef9e9575b4be843aa64478c4b1e0422a751`; the COSMIC v3.3.1 GRCh37 SBS table SHA256 is `9d9c212cc14aa660da1b7fd9e60e696c8e2c885a7606c27d20ad8997954554ea`.

## Full-campaign checks

| Dataset / calculation | Scope | Reference check |
| --- | --- | --- |
| Human COSMIC | 73,214 transcripts × 79 signatures × 1,000 paths = 5,783,906,000 | Complete checkpoint audit and checksummed export; zero incomplete paths |
| Two bat COSMIC controls | 31,185 Molossus + 23,213 Myotis transcripts × 8 signatures × 1,000 paths = 435,184,000 | Both cohorts audited and exported; zero incomplete paths |
| Six experimental profiles | All three cohorts × 6 profiles × 1,000 paths = 765,672,000 | Full audit and 18 dataset/profile summaries |
| Sampler validation | 79 × 10,000 sequences of 2,200 bp | 790,000 substitutions; all conditional-oracle checks passed |
| S11 historical comparison | 3,740 matched transcripts | Archived estimates and sequence identities verified |
| S12 graph | 1,000 complete paths, 3,152 states, 3,623 edges | Sequence reconstruction, chronological edges and stopping criterion reaudited |
| S12 breast-SNV preparation | 354 breast observations; 235 mapped observations; 174 distinct substitutions | Archive hashes, phenotype join, target coordinates and reference alleles checked |
| S12 graph annotation | 19 single-substitution allele matches | Numeric audit passed; supplied screenshot reviewed against the declared rule |

The human reference includes 19,309 unique genes and retains 73 internal-stop transcripts under the original selection policy. Bat counts are transcripts, not genes or one-to-one orthologues. Internal-stop counts and input selection policies are recorded by each preparation manifest.

## Numerical comparisons

The principal manuscript PCA fits the 60 signatures without sequencing artifacts. It mean-centres the matrix without variance scaling, uses full SVD and preserves the historical sign convention. Figure 4 displays the 41 known-aetiology signatures from that fitted model.

| Quantity | Preserved historical matrix | Checked full calculation |
| --- | ---: | ---: |
| PC1 + PC2 variance | 73.24703% | 73.23996% |
| PC1–PC10 variance | 90.20846% | 90.20196% |
| Transition fraction versus mean RUMS, R² | 0.71495226 | 0.71485790 |

The rounded 73%, 90% and R²=0.715 are retained. Full PC1–PC10 loading-vector cosines are at least 0.99923324. The displayed PC5 maximum changes from SBS10a to SBS10b, with a small coefficient difference (0.26022800 versus 0.26070322); both have similarly large loadings. PC8/SBS84, PC9/SBS22 and PC10/SBS1 remain the displayed maxima.

SBS2 has the highest human mean and across-transcript SD in each recorded signature group. Mean sampler cosine to the original profiles is 0.9751868632, compared with the historical manuscript value 0.974; available sequence opportunities affect realized channel frequencies, so an independent conditional oracle supplies the software check.

## Native-process species comparison

| Dataset / process | Mean RUMS | Median RUMS |
| --- | ---: | ---: |
| Human / COSMIC SBS2 | 0.789106 | 0.679 |
| Molossus / its APOBEC3A profile | 0.665150 | 0.648 |
| Myotis / its APOBEC3A3 profile | 0.787936 | 0.727 |

Myotis has the largest median; the human mean is slightly larger than Myotis. The preserved rank-test directions do not test an arithmetic-mean ordering. The human–Myotis rank-biserial effect is small (−0.06860 in the human-minus-Myotis orientation) despite a small p-value and large transcript count. Transcript isoforms and repeated signature measurements are not independent biological replicates; these analyses do not establish adaptation or replace orthologue-matched tests. Floating-point p-values of zero indicate numerical underflow.

## Supplementary scope

S11 uses archived estimates on the same 3,740 transcript identities for both criteria: mean 0.8116874332 under protein identity and 3.8720010695 under the historical BLOSUM criterion. It is historical postprocessing; a new corrected-package full-cohort BLOSUM sensitivity run is a separate experiment.

S12 uses a new 1,000-path SBS2/BLOSUM80 simulation (mean pre-terminal neutral steps 2.623, median 2.0). COSMIC v101/GRCh38 observations were mapped from the source transcript to ENST00000357654.9 using GENCODE v40 genomic intervals and reference-allele checks. Nineteen single-substitution nodes match reported breast SNVs. Another 288 multi-substitution nodes contain at least one reported component allele; this does not establish co-occurrence or observed full genotypes. The historical red-point rule is unrecovered, and the current caption explicitly defines the reconstructed panel.

BLOSUM80 distance is an alignment-score proxy, not experimentally measured functional loss. Exact historical random streams and the original experimental-profile probability CSVs are unavailable. The six profiles are reconstructed from the included motif-average workbook with a documented normalization policy. Wet-lab reads/upstream processing are outside these code workflows.

## Validation records and future runs

Local automated tests and wheel-installation checks verify software behavior. Full datasets are validated by run manifests, complete checkpoint audits, export hashes and numerical reports. GitHub Actions defines bounded software checks; consult its actual status for a particular commit before citing remote CI as passed. Preserve the exact exported figures and captions with the manuscript, including publication-scale typography checks. Reproducibility means a declared code/input/criterion configuration, not a claim of pixel-identical historical artwork or unchanged historical random streams.
