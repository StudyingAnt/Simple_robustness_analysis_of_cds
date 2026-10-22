# Recovered historical supplementary inputs

Additional archive files were located on 2026-10-21, after the initial reconstruction examples were prepared. This updates S11 input recovery and clarifies the surviving S12 source. The core validation programs are unchanged.

## S11: numerical source and transcript identities recovered

`Result_analysis.ipynb` reads 3,740 values from the `lrho` dataset of `lrho_sbs2.h5`, associates them with the first 3,740 records of `gencode.v40.pc_transcripts.nopary.cdsplus.longest.fa`, and reads protein-identity SBS2 values from `all_gencode_rums_profile.csv`. `run_l_robustness.py` specifies SBS2, the BLOSUM criterion, threshold 20 and 1,000 simulations per transcript; it does not record historical random seeds.

All 3,740 associated transcript IDs and sequences match the current pinned GENCODE v40 CDS-plus input. The archived `sbs2_partial.png` visually matches the manuscript S11 panel. The saved SVG contains 73,214 scatter observations under protein identity and 3,740 under BLOSUM80. Thus its two plotted populations differ. Matching-subset selection code also survives in the notebook, but the saved SVG did not plot equal population sizes.

The recovered BLOSUM results can be compared with identity results for the same 3,740 transcript names without new simulation. Their historical means are 3.8720010695 (BLOSUM) and 0.8116874332 (matched identity); every recovered BLOSUM estimate is at least its corresponding identity estimate. This checks archived estimates, not paired random paths or a new run under the corrected package.

The separate `robustness_sbs2_blosum_20.h5` contains 985 explicitly named CGC-selected results. It is a different exploratory dataset and must not replace the 3,740-result source above. The previously prepared twelve-transcript reconstruction is also a new, bounded example rather than the original S11 population. Prioritize the recovered sources for original-panel validation.

Use [the portable historical S11 analysis and SeaWulf wrapper](S11_HISTORICAL.md) to reproduce the matched-population postprocessing. The derived `source_data/s11_historical_values.csv.gz` includes both estimates and their transcript identities; the accompanying manifest fingerprints the archived HDF5, FASTA and notebook. Full original HDF5/FASTA files are not bundled and HDF5 software is not required by this analysis.

## S12: graph code and artwork recovered; observed-variant mapping unresolved

`Result_analysis.ipynb` selects BRCA1-203 and calls SBS2/BLOSUM threshold-20 simulation with `n_sim=1000`. It constructs a directed graph by deduplicating DNA sequence states, stores mutation annotations and scores on vertices, and saves a Kamada–Kawai graph. Additional exploratory layouts are present.

The recovered `brca-1_sbs2.svg` and PNG have no red observed-variant annotations, whereas the manuscript S12 image does. The inspected notebooks and checkpoint copies contain no breast-cancer variant import or observed-variant matching rule. `Cosmic_CancerGeneCensus_v99_GRCh38.tsv` is a gene census, not a BRCA1 mutation-observation table. The inspected `test.h5` contains one-step trajectories for ten other transcripts, not the BRCA1 graph.

The simulation graph source and historical artwork are therefore recovered. The table and rule connecting observed breast-cancer variants to red points remain unresolved within the inspected archive. A new COSMIC download requires a declared release, tissue filter, variant type, coordinate mapping and matching rule. It does not automatically recover the historical rule. The portable 100-path radial-layout example remains a separate reconstruction.

**Subsequent reconstruction:** COSMIC v101/GRCh38 archives were downloaded, verified and filtered to 174 distinct breast SNVs in the simulated BRCA1 target coordinates. The [full-data graph annotation](BRCA1_OBSERVED_SNVS.md) passed on HPC, with 19 single-substitution allele matches and 288 multi-substitution nodes containing a reported component allele. The author-supplied screenshot has also been visually reviewed alongside the historical manuscript image. The regenerated panel restricts red points to single-substitution nodes, whereas historical red points span several layers. Retain this changed matching definition in the caption; the historical source-table release and original red-point rule remain unrecovered.

## Fingerprints

The archive originals are unchanged; local research copies and a checksum manifest preserve the recovered files. Only the derived historical estimates and transcript IDs/names are bundled for this stage; the raw archive is not published.

| Source file | SHA256 |
|---|---|
| `Result_analysis.ipynb` | `abbc7b8545dad4fa05b0f4234bdae2ca8a8fa4336b37f2df0baa1151052c54a1` |
| `lrho_sbs2.h5` | `9777028bb2b855aa3ecf78adea6247e1203b7f9c8967cc1d09bd34af099abe8d` |
| `gencode.v40.pc_transcripts.nopary.cdsplus.longest.fa` | `6e576d0d49d47479587309bb45ef49fe8995fe80163a666a838a9cfc47eaa978` |
| `sbs2_partial.svg` | `6ecdaf1e70ac46f0b1cc1798b5bee7a2db5cd07c93192d4487afff6d74ca5844` |
| `brca-1_sbs2.svg` | `c9025956e0c7bb870883bcc9f11ba3df4a961e3890d095c0be9ee5911bfe628f` |
