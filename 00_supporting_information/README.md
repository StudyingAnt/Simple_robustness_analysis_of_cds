# Supporting tables and calculated RUMS data

Download the completed numerical matrices and current supplementary tables here
without rerunning Sinabro or using an institutional HPC account. These are the
audited calculations, not synthetic examples: 1,000 trajectories per
transcript/profile combination and master seed 42.

## Supplementary tables

| File | Contents |
| --- | --- |
| [Table S1](Table_S1.xlsx) | Human repeated-measures analysis: 19,309 fixed gene representatives, 79 profiles, 3,081 paired Wilcoxon contrasts/Holm correction, medians, paired median differences/pointwise 95% intervals and rank effects |
| [Table S2](Table_S2.xlsx) | SBS aetiology categories, with corrected internal title |
| [Table S3](Table_S3.xlsx) | Descriptive human/bat CDS-cohort distributions and contrasts under eight fixed COSMIC profiles |
| [Table S4](Table_S4.xlsx) | Descriptive cohort distributions and contrasts under six experimental APOBEC profiles and cohort-specific profiles |
| [Table S5](Table_S5.xlsx) | High-RUMS transcript outliers, gene annotations, CDS length, total GC and GC3 |
| [Table S6](Table_S6.xlsx) | All six experimental NNC preferences: absolute motif-average percentages and normalized per-opportunity weights |

Tables S3 and S4 were regenerated from the same audited exports. Their
`Shared_SBS2` sheets are identical. These revised tables replace the old
transcript-level species significance interpretation with descriptive cohort
medians, differences and rank effects. No species-level confidence intervals are
claimed: transcripts are not independent animals, and the bat inputs are not
verified one-to-one orthologue sets. Table S1 instead uses a paired human gene
design. Numeric p-value underflow is identified with log-scale/display columns;
a stored zero does not mean an exactly zero probability.

The current Tables S1–S6 are directly in this directory. Earlier source
workbooks are preserved in [historical_source_tables/](historical_source_tables/)
and do not use the current manuscript numbering.

## Download the calculated matrices

All are gzip-compressed CSVs with headers and explicit identifiers. Values are
mean accepted synonymous substitutions before the terminal nonsynonymous event,
including zero-count paths, not mutation rates or measured functional loss.
COSMIC profiles use the archived v3.3.1 GRCh37 table; bat controls contain eight
profiles. Experimental matrices apply the same six NNC profiles to every cohort.

| Matrix | Rows × profiles |
| --- | --- |
| [Human transcript × COSMIC](data_matrices/human_cosmic_rums_mean.csv.gz) | 73,214 × 79 |
| [Human gene representative × COSMIC](data_matrices/human_gene_representatives_cosmic_rums_mean.csv.gz) | 19,309 × 79 |
| [Molossus transcript × COSMIC](data_matrices/molMol2_cosmic_rums_mean.csv.gz) | 31,185 × 8 |
| [Myotis transcript × COSMIC](data_matrices/myoMyo6_cosmic_rums_mean.csv.gz) | 23,213 × 8 |
| [Human transcript × experimental APOBEC](data_matrices/human_apobec_rums_mean.csv.gz) | 73,214 × 6 |
| [Human gene representative × experimental APOBEC](data_matrices/human_gene_representatives_apobec_rums_mean.csv.gz) | 19,309 × 6 |
| [Molossus transcript × experimental APOBEC](data_matrices/molMol2_apobec_rums_mean.csv.gz) | 31,185 × 6 |
| [Myotis transcript × experimental APOBEC](data_matrices/myoMyo6_apobec_rums_mean.csv.gz) | 23,213 × 6 |

Human gene matrices retain one longest CDS per stable gene, with transcript
accession tie-breaking before examining RUMS. `gene_id` and
`representative_transcript` are separate columns. Selection records and all
three transcript annotation files are included. Bat rows remain **transcripts**;
their frozen `mRNA` identifiers must not be presented as verified orthologous
genes.

Plain Python can read the values after cloning this repository:

```python
import csv
import gzip

with gzip.open(
    "00_supporting_information/data_matrices/human_cosmic_rums_mean.csv.gz",
    "rt", newline=""
) as stream:
    reader = csv.DictReader(stream)
    print(next(reader))
```

## Six NNC preferences

[Table S6](Table_S6.xlsx) has a sheet per enzyme.
[profile_weights.csv](apobec_preferences/profile_weights.csv) gives
all 16 × 6 normalized weights together. `source_values.csv` and `profiles.json`
retain absolute percentages and exact simulation inputs. Weights sum to one per
profile and are applied per eligible site, followed by normalization over events
available in the current sequence; they are not whole-genome event probabilities.
Only C-to-T is represented. Motif occurrences do not supply independent reaction
replicate counts, background rates or replicate error bars. These are the
workbook-derived inputs used for the audited recalculation; the historical
probability CSV was not recovered.

## Provenance and checks

`data_manifest.json` records input/configuration hashes, source-export hashes,
dimensions and compressed/uncompressed file hashes. `SHA256SUMS` lists every
data/workbook checksum. The standard-library verifier checks hashes, dimensions,
identifiers and profile weights without generating any trajectories:

```bash
python 00_supporting_information/verify_data.py
```

These numerical files can accompany the paper as supplementary data as well as
being available on GitHub. A commit-specific link identifies the exact snapshot;
no DOI-bearing deposit is claimed. Full checkpoints and reference genomes are
not included. Restricted COSMIC patient observations and credentials are not
part of this bundle.
