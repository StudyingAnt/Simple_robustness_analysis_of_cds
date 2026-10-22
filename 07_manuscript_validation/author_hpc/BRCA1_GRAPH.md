# S12: fresh BRCA1 simulation graph

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

This stage executes 1,000 fresh SBS2/BLOSUM paths for the original BRCA1-203 case using Sinabro 3.0.0. It is separate from the completed historical S11 postprocessing. It requires the existing pinned human FASTA and package checkout. Actual COSMIC mutation observations can be acquired independently and annotated later; they are not required for this simulation.

## Completed full HPC verification

The author supplied a successful full HPC report on 2026-10-21: **1,000 complete paths, zero incomplete paths**, graph audit passed, **3,152 DNA states** and **3,623 single-substitution edges**. Mean neutral steps are **2.623**, median **2.0**. Reported elapsed time was **11.904 seconds**. The four scientific data files (path counts, graph nodes, chronological edges and reference FASTA) have exactly the same reported SHA256 values as the verified local run. The supplied report did not include its Slurm job ID or absolute run directory.

The initial HPC simulation-only image hashes differ from the local hashes; their cause was not established. The subsequent COSMIC annotation passed, and the author-supplied screenshot of the regenerated panel has been reviewed as described in [the annotation guide](../BRCA1_OBSERVED_SNVS.md). This does not recover the historical red-point rule or establish byte-identical artwork. Existing successful simulations do not need resubmission for annotation or documentation updates.

## Submit one job on SeaWulf

For the independent observed-variant input, see [COSMIC v101/GRCh38 download instructions](../../01_data_downloading_and_preprocessing/COSMIC_VARIANTS.md). The portable downloader checks official archive SHA256 values and prompts for account credentials only in an interactive terminal. v101 is a date-based candidate, not a confirmed historical release. This acquisition does not require another simulation job.

**Observed-input update:** the author-reported HPC extraction passed for 1,984,965 source rows and retained 3,037 BRCA1 observations, all annotated on **ENST00000471181.7**. Local GENCODE v40 sequence checks show this source CDS is 5,655 nt, whereas the simulated **ENST00000357654.9** CDS is 5,592 nt. Source CDS coordinates cannot be reused directly. The completed preparation therefore joins breast phenotypes and uses GRCh38 genomic-to-target-CDS mapping with reference-allele checks; the existing graph needs no new trajectories.

The completed [breast selection and coordinate-mapping stage](../../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md) used the two raw TARs and a frozen 22-segment GENCODE v40 map whose complete 5,592-nt CDS was independently verified against NC_000017.11 genomic sequence.

**Mapping and annotation completed:** the author supplied passing full reports: 354 breast observations, 235 accepted mapped observations, **174 distinct target SNVs**, no attention reasons. Both official archive hashes, FASTA, map and code hashes matched. The [graph allele-matching stage](../BRCA1_OBSERVED_SNVS.md) passed using the existing 1,000 paths: **19 single-substitution allele-match nodes**, **288 multi-substitution nodes containing a reported component allele**. Component matches do not establish observed multi-mutant genotypes. Author-screenshot content/layout review is complete; the caption must retain the changed matching definition.

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs
sbatch Simple_robustness_analysis_of_cds/07_manuscript_validation/author_hpc/run_brca1_graph_seawulf.slurm
```

The job requests eight shared CPUs, 8 GB and a two-hour limit and activates the existing `PROJECT/sinabro-revision` environment. The requested limit is not an execution-time estimate. No new dependency installation is needed. Only BRCA1 is simulated, not the full human or bat datasets.

The dedicated wrapper uses **1,000** paths, matching the count in the recovered original notebook. The older generic supplementary wrapper retains a separate 100-path S12 example; use the dedicated wrapper above for this full case.

## Scientific configuration

- BRCA1-203, **ENST00000357654.9**, 5,592 coding nucleotides plus the original one-base flanks.
- COSMIC v3.3.1, GRCh37 SBS2 profile, explicit signature checksum recorded. The signature profile's reference label does not change the pinned GENCODE v40 transcript sequence.
- BLOSUM80 global alignment, gap-open −10, gap-extension −0.5; stop at alignment-score distance **≥20**, `stop_policy=score`.
- Master seed42, independent deterministic transcript/path-keyed streams; the historical exact streams were not recovered.
- Safety cap10,000 proposed mutations per path. Any incomplete path fails the run rather than contributing a truncated result.

## Outputs and validation

Results are intermediate simulation records under:

```text
PROJECT/runs/brca1_sbs2_blosum20_n1000_JOBID/
```

- **`supplementary_summary.json`**: scientific configuration, complete-path/graph audit status, counts, reference/profile/source hashes and descriptive path summaries.
- `trajectory_counts.csv`: path index, seed, neutral-step count and terminal reason/distance.
- `s12_reference.cdsplus.fa`: exact reference transcript with CDS bounds.
- `s12_graph_nodes.csv`: unique sequence-state hashes, reference-relative CDS substitutions, Hamming distance, BLOSUM distance and diagnostic plot coordinates.
- `s12_graph_edges.csv`: source/target nodes, path and chronological mutation step. This supports reconstructing each path even when states are shared or revisited.
- `figureS12_simulation_only.pdf/.png`: diagnostic plot; not the final observed-variant panel.

The audit independently reconstructs every DNA state from the reference and its CDS substitutions and verifies its hash. It checks each edge is exactly one nucleotide change, each path starts at the original reference, chronological steps and accepted mutation counts agree, and all paths stop at the declared score threshold. All nodes must belong to the audited paths.

After completion, share `supplementary_summary.json`. Expected completion fields:

```json
{
  "simulation_passed": true,
  "fresh_simulation": true,
  "corrected_package_blosum_recalculation": true,
  "n_trajectories": 1000,
  "graph_audit_passed": true,
  "complete_paths": 1000,
  "incomplete_paths": 0
}
```

`original_figure_reproduced` and `observed_cosmic_variant_overlay_reproduced` remain **false** for this intermediate stage. Successful simulation is separate from recovering and checking the actual observed variants. The generated DNA records can be annotated later without new trajectories. The original observation-to-node rule is unresolved, so any replacement annotation must state its own evidence and matching definition rather than claiming the historical rule was recovered.

The diagnostic layout uses deterministic angles and a radius equal to actual reference DNA Hamming distance. The recovered original notebook used a Kamada–Kawai graph and also explored other layouts. These layouts are not claimed identical. BLOSUM distance is an alignment proxy, not experimentally measured functional loss.

## Portable execution

From the analysis repository root with the pinned environment active:

```bash
python 07_manuscript_validation/run_supplementary_cases.py \
  --stage s12 \
  --input /PATH/TO/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --sinabro-source /PATH/TO/sinabro \
  --out-dir /PATH/TO/PROJECT/runs/brca1_sbs2_blosum20_n1000 \
  --n-traj 1000 --workers 8 --seed 42 --max-steps 10000
```

All output directories must be new. The Slurm defaults can be overridden with `PROJECT`, `ANALYSIS_DIR`, `SINABRO_ENV`, `SINABRO_SOURCE`, `SINABRO_FASTA` and `SINABRO_BRCA1_OUTPUT`.
