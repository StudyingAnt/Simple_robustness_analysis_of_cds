# Six experimental APOBEC profiles: simulation before plotting

This workflow applies each of the six measured bat enzyme preferences to the same human, Molossus molossus and Myotis myotis CDS cohorts used by the revised COSMIC calculations. These are empirical NNC preferences, distinct from COSMIC SBS2/SBS13. Python programs are portable; only the Slurm wrappers contain SeaWulf settings.

| Dataset | Transcripts | Enzyme profiles | Paths per combination | Total paths |
| --- | ---: | ---: | ---: | ---: |
| Human GENCODE v40 | 73,214 | 6 | 1,000 | 439,284,000 |
| HLmolMol2 | 31,185 | 6 | 1,000 | 187,110,000 |
| HLmyoMyo6 | 23,213 | 6 | 1,000 | 139,278,000 |
| Total | 127,612 | 6 | 1,000 | 765,672,000 |

The full experimental campaign has not been run by the assistant. Existing COSMIC outputs are reused as separate controls and are not recomputed here. Plotting and comparative statistics follow the full experimental archive audit.

## Source workbook and reconstruction policy

The original private workbook is `C to T Mutation data of six bat enzymes.xlsx`, SHA256 `af2e8dad39b7ebb3e36d4e55d388c7550967610398b7361dc06cdddaae5add5a`. Its values and generated probability tables are excluded from Git. The supplied Drive link required sign-in in the download check; transfer the original file privately rather than relying on an anonymous HPC download.

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

## SeaWulf: run one stage at a time

Use the existing revision environment (Python 3.9.15, NumPy 1.21.6, Biopython 1.79), existing extracted CDS inputs, and the refactored Sinabro package source. All sbatch commands below are issued from the HPC project root so relative log paths resolve correctly.

### 1. Update code and copy the workbook once

On HPC:

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs data/raw/bat_apobec
```

Transfer the original XLSX using the same file-transfer connection used for this HPC account. The required destination is:

```text
/gpfs/projects/DamaghiGroup/projects/sinabro_pub/data/raw/bat_apobec/C to T Mutation data of six bat enzymes.xlsx
```

For example, from a **local WSL terminal**, replacing `YOUR_HPC_LOGIN_HOST` with your usual SSH login hostname:

```bash
scp '/home/ceel/projects/sinabro_all/Simple_robustness_analysis_of_cds/data/raw/bat_apobec/C to T Mutation data of six bat enzymes.xlsx' \
  joonsong@YOUR_HPC_LOGIN_HOST:/gpfs/projects/DamaghiGroup/projects/sinabro_pub/data/raw/bat_apobec/
```

The source path is local to the assistant's current WSL workspace. Other computers should use their own copy of this workbook.

### 2. Prepare and verify profiles and inputs

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/experimental_apobec/prepare_experimental_apobec_seawulf.slurm
```

Check `logs/apobec-prepare-JOBID.out` and `.err` before proceeding. Preparation checks the frozen workbook, sample identities, motifs, percentages, human FASTA hash, bat extraction manifests and population counts. The final JSON should show `preparation_passed: true`, `combinations: 765672` and `production_paths_at_n1000: 765672000`.

Default generated directories:

```text
data/processed/experimental_apobec_v1/
  profiles.json
  source_values.csv
  preparation_manifest.json
data/processed/experimental_apobec_inputs_v1/
  inputs_manifest.json
```

Preparation refuses to overwrite an existing campaign input directory. Existing prepared profiles are verified before use. For a separate preparation, set `SINABRO_APOBEC_PROFILES` and `SINABRO_APOBEC_INPUTS` to new locations and keep those environment variables consistent in later stages. `SINABRO_APOBEC_WORKBOOK`, `SINABRO_HUMAN_FASTA` and `SINABRO_BAT_PROCESSED` override input locations. Raw and processed data remain outside both Git checkouts.

### 3. Check execution on HPC

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/experimental_apobec/validate_experimental_pilot_seawulf.slurm
```

This tests 12 length-stratified transcripts per dataset, all six profiles and 100 paths per combination, serially and with eight workers: 43,200 simulated paths in total. Each mode is independently audited and every path count/termination code is compared. The pilot omits internal-stop transcripts only from this bounded check; the full campaign retains them.

Require `validation_passed: true` in:

```text
runs/experimental_apobec_pilot_JOBID/experimental_pilot_report.json
```

### 4. Submit the complete campaign once

Replace `PILOT_JOBID` with the passing pilot's actual job ID:

```bash
export SINABRO_APOBEC_PILOT="$PWD/runs/experimental_apobec_pilot_PILOT_JOBID"
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/experimental_apobec/run_experimental_apobec_seawulf.slurm
```

One submission contains **50 array tasks**, with at most eight tasks running concurrently. Each task requests one node and 96 workers; the maximum concurrent request is eight nodes and 768 workers. Allocation depends on scheduler availability. Every task runs six dataset/shard cases sequentially: shard k on all three datasets, then shard k+50 on all three datasets. Together these cover 100 shards per dataset exactly once, with no later range to submit manually. The four-hour limit applies to each allocation and is a requested limit, not a measured completion estimate.

Production requires the passing HPC pilot with matching input/profile hashes, simulation source, profile-loader source, package source, master seed and step cap. The full and pilot path counts intentionally differ. Keep code and inputs fixed during the campaign.

Default output:

```text
runs/hpc_experimental_apobec_n1000_seed42/
  human/{run_manifest.json,transcripts/,shards/}
  molMol2/{run_manifest.json,transcripts/,shards/}
  myoMyo6/{run_manifest.json,transcripts/,shards/}
```

Each transcript archive contains all six profiles and every path's count and termination code. Seeds include master seed, dataset, transcript ID and enzyme ID, so scheduling and worker counts do not alter path results. Checkpoints are published atomically and validated before resume. A duplicate live shard is rejected by a lock. Configuration mismatches require a new output directory; existing results are not overwritten.

For an interrupted task, resubmit only its original task index using `sbatch --array=INDEX .../run_experimental_apobec_seawulf.slurm` with the same inputs, pilot and run settings. Completed transcript checkpoints in that task are checked and skipped. Do not change the task-count or shard plan. A stored incomplete path is not silently treated as a successful result or repaired on resume; investigate it and use a separate, explicitly changed run if necessary. Existing unrelated queued jobs may still affect the account's submission limit even with this bounded array.

### 5. Audit the completed archive

After all production array tasks finish:

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/experimental_apobec/validate_experimental_results_seawulf.slurm
```

Require `stage_passed: true` in `runs/experimental_apobec_audit_JOBID/experimental_audit_manifest.json`. The audit reads every archive and checks populations, configuration fingerprints, all shard completion records, individual path completion and stored statistics. It performs no new simulation. Table export, plots and statistical comparisons are subsequent stages.

## Other computers

Activate the revision environment and use the portable programs directly. Paths below are examples; supply the local copies of the same frozen inputs:

```bash
python 06_Bat_genome/experimental_apobec/prepare_profiles.py \
  --workbook '/PATH/TO/PROJECT/data/raw/bat_apobec/C to T Mutation data of six bat enzymes.xlsx' \
  --out-dir /PATH/TO/PROJECT/data/processed/experimental_apobec_v1
python 06_Bat_genome/experimental_apobec/prepare_inputs.py \
  --human-input /PATH/TO/PROJECT/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa \
  --bat-root /PATH/TO/PROJECT/data/processed/bat1k_historical_v1 \
  --profiles-dir /PATH/TO/PROJECT/data/processed/experimental_apobec_v1 \
  --out-dir /PATH/TO/PROJECT/data/processed/experimental_apobec_inputs_v1
python 06_Bat_genome/experimental_apobec/validate_experimental_pilot.py \
  --inputs-dir /PATH/TO/PROJECT/data/processed/experimental_apobec_inputs_v1 \
  --out-dir /PATH/TO/PROJECT/runs/experimental_apobec_pilot \
  --sinabro-source /PATH/TO/PROJECT/sinabro --workers 4 --n-traj 100
```

`run_experimental_apobec.py --help` provides a single-dataset runner. `run_campaign.py --help` provides the allocation-independent task driver; call all task indices for the chosen plan. Input manifests contain local absolute paths, so prepare them on each machine rather than copying a manifest with paths from another machine.

## Local validation before handoff

The original workbook and all three actual cohorts passed preparation locally on 2026-10-19. The complete bounded pilot passed: 43,200 paths, all six independent archive audits passed, serial/four-worker results identical path by path, four simultaneous cases observed and no incomplete paths. An independent sampler check drew 2,000 events per profile and checked event eligibility, substitutions and observed counts against separately computed per-position probabilities. A further bounded production-driver check completed six dataset/shard cases (132 transcripts, 1,584 paths), including the longest input per dataset, and independent checkpoint auditing passed; resume skipped all completed transcripts and preserved every checkpoint hash. The analysis repository's 104 unit tests passed, including source corruption, identities, normalization, immutable outputs, checkpoint statistics, incomplete paths, pilot provenance and the complete 50-task allocation plan. All four new Slurm files passed bash syntax checks. These checks establish bounded implementation validation, not full experimental reproduction.
