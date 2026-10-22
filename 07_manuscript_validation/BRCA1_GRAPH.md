# S12: BRCA1 mutation graph

Run from the analysis repository root with the pinned environment and package checkout available:

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
SINABRO_SOURCE=/ABSOLUTE/PATH/TO/YOUR/sinabro

python 07_manuscript_validation/run_supplementary_cases.py \
  --stage s12 \
  --input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --sinabro-source "$SINABRO_SOURCE" \
  --out-dir "$PROJECT/runs/brca1_sbs2_blosum20_n1000" \
  --n-traj 1000 --workers 1 --seed 42 --max-steps 10000
```

This simulates only BRCA1-203 (**ENST00000357654.9**, 5,592 coding nt) under SBS2, using BLOSUM80 global alignment with gap-open −10, gap-extension −0.5 and stopping distance **≥20** (`stop_policy=score`). The fixed COSMIC v3.3.1 GRCh37 profile label does not change the GENCODE v40 reference sequence.

The new output directory contains `supplementary_summary.json`, `trajectory_counts.csv`, `s12_reference.cdsplus.fa`, `s12_graph_nodes.csv`, `s12_graph_edges.csv` and diagnostic `figureS12_simulation_only.pdf/.png`. The graph audit reconstructs DNA states and checks every edge is a single substitution, path chronology and termination. A passing full run has 1,000 complete paths, no incomplete paths and `graph_audit_passed=true`.

For the checked source/environment/seed, the reference graph has **3,152 states and 3,623 edges**. Diagnostic radius is reference DNA Hamming distance; angles are deterministic. This layout and random streams are not claimed identical to historical artwork. BLOSUM score distance is an alignment proxy rather than measured functional loss.

Next, [prepare reported breast SNVs](../01_data_downloading_and_preprocessing/COSMIC_BREAST_VARIANTS.md) and [annotate the existing graph](BRCA1_OBSERVED_SNVS.md). Annotation needs no new simulation. Historical COSMIC release and red-node semantics were not recovered; the regenerated panel states its current matching rule.

The authors ran this calculation on SeaWulf. Their [actual wrapper, configuration and execution record](author_hpc/BRCA1_GRAPH.md) are separate reference material. SeaWulf requires institutional access; the Python command above does not. Adapt the optional wrapper to your own site if you use HPC.
