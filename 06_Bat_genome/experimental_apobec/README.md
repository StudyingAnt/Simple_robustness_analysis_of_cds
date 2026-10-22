# Six experimental APOBEC profiles

Apply the six measured bat enzyme preferences to the human, Molossus and Myotis CDS cohorts. These empirical NNC preferences are distinct from COSMIC SBS2/SBS13. Use the [pinned environment](../../CONDA_SETUP.md), [prepared human input](../../01_data_downloading_and_preprocessing/README.md) and [both prepared bat inputs](../README.md).

Run the following commands from the analysis repository root. They use ordinary Python, with no scheduler, institutional account or pilot report required.

## 1. Prepare the profiles and input manifest

```bash
PROJECT=/ABSOLUTE/PATH/TO/YOUR/sinabro_project
WORKERS=1

python 06_Bat_genome/experimental_apobec/prepare_profiles.py \
  --workbook "06_Bat_genome/experimental_apobec/source_data/C to T Mutation data of six bat enzymes.xlsx" \
  --out-dir "$PROJECT/data/processed/experimental_apobec_v1"

python 06_Bat_genome/experimental_apobec/prepare_inputs.py \
  --human-input "$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa" \
  --bat-root "$PROJECT/data/processed/bat1k_historical_v1" \
  --profiles-dir "$PROJECT/data/processed/experimental_apobec_v1" \
  --out-dir "$PROJECT/data/processed/experimental_apobec_inputs_v1"
```

The workbook is included in Git. Profile preparation writes `profiles.json`, `source_values.csv` and `preparation_manifest.json`. Input preparation checks the frozen source sequences and produces `inputs_manifest.json`. Each output directory must be new.

## 2. Run all six profiles on each dataset

```bash
HUMAN_INPUT="$PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa"
BATS="$PROJECT/data/processed/bat1k_historical_v1"

for DATASET in human molMol2 myoMyo6; do
  if [ "$DATASET" = human ]; then
    INPUT="$HUMAN_INPUT"
  else
    INPUT="$BATS/$DATASET/$DATASET.cdsplus.fa"
  fi
  python 06_Bat_genome/experimental_apobec/run_experimental_apobec.py \
    --input "$INPUT" --dataset "$DATASET" \
    --profiles-dir "$PROJECT/data/processed/experimental_apobec_v1" \
    --out-dir "$PROJECT/runs/experimental_apobec_n1000_seed42/$DATASET" \
    --n-traj 1000 --seed 42 --max-steps 10000 --workers "$WORKERS"
done
```

The full campaign contains **765,672,000 paths**: 439,284,000 human, 187,110,000 Molossus and 139,278,000 Myotis. It is a large manuscript calculation, not a short workstation demonstration. Choose the worker count to suit your computer. The direct runner needs the prepared profiles but does not require the authors' serial/parallel pilot or task-allocation controller.

## 3. Audit and export the existing results

```bash
python 06_Bat_genome/experimental_apobec/validate_campaign.py \
  --inputs-dir "$PROJECT/data/processed/experimental_apobec_inputs_v1" \
  --run-root "$PROJECT/runs/experimental_apobec_n1000_seed42" \
  --out-root "$PROJECT/runs/experimental_apobec_audit" \
  --n-traj 1000 --workers "$WORKERS"

python 06_Bat_genome/experimental_apobec/export_results.py \
  --inputs-dir "$PROJECT/data/processed/experimental_apobec_inputs_v1" \
  --run-root "$PROJECT/runs/experimental_apobec_n1000_seed42" \
  --audit-root "$PROJECT/runs/experimental_apobec_audit" \
  --out-root "$PROJECT/runs/experimental_apobec_export" \
  --n-traj 1000 --workers "$WORKERS"
```

Inspect `experimental_audit_manifest.json` and `export_manifest.json`; both should report `stage_passed=true`. The export includes `human/`, `molMol2/`, `myoMyo6/` and `profile_comparison_summary.csv` (18 dataset/profile rows). Continue with [manuscript statistics and panels](../../07_manuscript_validation/README.md). Output directories must be new; full audit/export uses the full frozen populations.

## Source workbook and reconstruction policy

The original workbook is versioned at [source_data/C to T Mutation data of six bat enzymes.xlsx](source_data/C%20to%20T%20Mutation%20data%20of%20six%20bat%20enzymes.xlsx), SHA256 `af2e8dad39b7ebb3e36d4e55d388c7550967610398b7361dc06cdddaae5add5a`. The project owner authorized publication of this source file for reproduction. Git pull retrieves the workbook together with the code; no Drive sign-in or separate file transfer is required. Generated probability tables and analysis outputs remain excluded from Git.

| Source sample | Profile ID |
| --- | --- |
| B7 | HLmolMol2_APOBEC3A |
| B8 | HLmolMol2_APOBEC3C |
| B24 | HLmolMol2_APOBEC1 |
| B22 | HLmyoMyo6_APOBEC3A3 |
| B23 | HLmyoMyo6_APOBEC3C |
| B27 | HLmyoMyo6_APOBEC1 |

The Myotis source label APOBEC3A3 is preserved. Historical aggregate CSVs called that column APOBEC3A.

The workbook identifies the 16 NNC values as motif averages. For each enzyme, the per-opportunity weight is the motif's mean percentage divided by the sum of its 16 mean percentages. The parser verifies that percentage equals 100 times mean C-to-T count divided by the enzyme's sequenced-molecule total. Normalizing those mean counts gives the same weights. Substrate motif occurrence counts are recorded for provenance and are **not multiplied into these weights**.

The historical package's `mut_types` sampler assigned each eligible sequence occurrence its supplied motif weight. This reconstruction follows that interpretation. The original `molMol2_apobec.csv` and `myoMyo6_apobec.csv` conversion tables have not been recovered, so exact identity to the historical probability tables is not established. Original aggregate RUMS tables are not a substitute for those inputs. The preparation manifest records this limitation explicitly.

The model uses `NN[C>T]` with the reverse-complement G-to-A events, `strand=both` and equal strand factors of 0.5. Event probabilities are normalized over all available weighted occurrences; equal strand factors do not force half the sampled events onto each strand. C targets require two upstream bases and reverse-complement G targets require two downstream bases. The existing one-base CDS flanks are retained, targets stay inside the CDS, and a target without its full context is excluded. No extra flanks or centered SBS96 contexts are invented.

This models relative sequence preference, not absolute enzyme activity or biological elapsed time. No uncertainty estimate for experimental replicates is inferred from motif means. Paths terminate at the first amino-acid sequence change. Saved counts record accepted neutral steps before that event. The 10,000-step cap and unavailable-event terminations are explicit incomplete paths; the full audit must pass before interpretation. The existing internal-stop filtering policy is retained for the full cohorts (73 human, 383 Molossus, 79 Myotis).


## Author execution

The authors performed the full campaign on Stony Brook University's SeaWulf HPC system, which requires institutional access. Their actual Slurm wrappers, allocation commands, pilot checks and completed-run information are retained under [author_hpc/](author_hpc/README.md). Researchers using an HPC system should adapt those examples to their own site's scheduler and environment. The public Python workflow above is independent of SeaWulf.
