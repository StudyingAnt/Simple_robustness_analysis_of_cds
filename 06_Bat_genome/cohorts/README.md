# Frozen historical transcript cohorts

These identifier lists come from the author's historical molMol2_rums_profile.csv
and myoMyo6_rums_profile.csv result tables, not from a new biological selection.
manifest.json records their source-table hashes, identifier-list hashes and counts:
31,185 Molossus molossus and 23,213 Myotis myotis transcripts. Sorting is lexical.
Only transcript IDs are distributed here; experimental enzyme data and CGC data
are not included.

Every listed ID matches the corresponding official allAnnotatedIsoforms BED
snapshot in bat1k_pilot_sources.json. For every listed record, chromStart equals
thickStart and chromEnd equals thickEnd. The extractor enforces this condition,
so clipping exon blocks to coding coordinates preserves these historical intervals.
It does not add UTR-bearing transcripts that the historical code may have excluded.

Reapplying the original all-exon start/stop/length filter to the current references
accepts exactly the Molossus list, but two additional Myotis transcripts:
mRNA9113 and mRNA9159. Both contain N bases in their CDS. Neither occurs in the
historical result table, and neither is added to the frozen reproduction cohort.
The exact historical reason for their absence is not recorded. For every listed
transcript, the new sequence agrees with that legacy extraction after uppercasing
the flanks. The Molossus legacy script left flanks in mixed case for 2,638 accepted
records; the new extractor normalizes all bases to uppercase.

Counts are transcripts, not genes. There is no representative-isoform selection or
one-to-one orthologue restriction. IDs can recur between species and must always
be interpreted within their species. Exact historical input FASTAs have not been
recovered; matching cohorts and new sequence validation do not establish that an
old simulation can be reproduced bit for bit.
