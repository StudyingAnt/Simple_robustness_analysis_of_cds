# Bat APOBEC motif-mean figure

This is a presentation update for manuscript Figure 6. It uses the existing
six-enzyme experimental input values without changing any RUMS calculation.
The older collaborator Methods document calls the same figure Figure 7;
the current manuscript identifies it as Figure 6.

Run with the repository's existing pinned environment:

```bash
python Extra_analysis/02_bat_assay_figure/redraw_assay_profiles.py \
  --profiles-dir /PATH/TO/prepared_experimental_apobec_profiles \
  --out-dir /PATH/TO/new_figure6_output
```

The input directory is the output of
`06_Bat_genome/experimental_apobec/prepare_profiles.py`: `profiles.json`,
`source_values.csv` and `preparation_manifest.json`. Both prepared data files
are authenticated against the manifest. Outputs are a PDF, PNG, caption draft
and provenance report. The directory must be new and outside the input data.
No scheduler or institutional HPC access is required.

All six panels use the same percentage axis. Bars retain the original absolute
motif-average percentages rather than the relative weights normalized across
16 motifs for simulation. NNC denotes two upstream bases and the terminal
cytosine that undergoes the C-to-T change.

The supplied assay description defines averages over occurrences of a motif
within the 369-bp DNA substrate. Those sequence positions are not documented
as independent enzyme reactions. The supplied records do not specify independent
reaction replicate n, replicate-level variability or a negative-control background
rate. Consequently, the script does not fabricate replicate error bars or assume
that missing replicate metadata means n=1. The caption is a draft and explicitly
records the missing information for author review.

The following would resolve the remaining Figure 6 reporting questions:

1. How many independent enzyme reactions or independently prepared libraries
   contributed to each sample, and are the per-repeat values available?
2. Was an enzyme-free or inactive-enzyme control sequenced, and were the reported
   percentages background-subtracted?

Read counts, motif-copy counts and the 16 different motif categories must not be
substituted for the number of independent assay repeats. A common-axis redraw
addresses the axis presentation request, not the remaining replication/control
information request.
