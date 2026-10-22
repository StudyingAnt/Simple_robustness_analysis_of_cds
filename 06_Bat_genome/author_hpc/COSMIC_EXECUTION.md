# Full-cohort bat COSMIC execution, audit and table export

> Author execution record: these commands and paths describe the authors’ SeaWulf setup. SeaWulf requires institutional access. For ordinary terminal/Python commands, use the stage README one level above. Researchers using an HPC system must adapt these examples to their own site.

Continue here after the CDS extraction and serial/parallel pilot in
HPC_REPRODUCTION.md have passed. The package and existing human runner/exporter
are unchanged. Python scripts are portable; only Slurm wrappers target SeaWulf.

## Scope and resources

| Setting | Value |
| --- | --- |
| Molossus transcripts | 31,185 |
| Myotis transcripts | 23,213 |
| COSMIC profiles, in historical order | SBS1, SBS2, SBS3, SBS4, SBS6, SBS9, SBS13, SBS18 |
| Trajectories per transcript/profile | 1,000 |
| Total trajectories | 435,184,000 |
| Master seed | 42 |
| Maximum steps per trajectory | 10,000; capped paths fail the final audit |
| Scientific shards | 100 per species, 200 species/shard cases |
| Submitted array jobs | 50; one submission covers the entire analysis |
| Work per array job | Two shard numbers, each processed for both species sequentially |
| Maximum simultaneous array elements | 8 |
| Resources per element | One node, 96 workers, four-hour wall limit |
| Maximum concurrent allocated workers | 768 across up to 8 nodes |

The four-hour limit is per array element, not a whole-run time estimate. Actual
allocation depends on the scheduler. Numerical threads are limited to one per
worker. Allocation task 0 processes Molossus shard 0, Myotis shard 0, Molossus
shard 50 and Myotis shard 50, in that order. Task 1 handles shard numbers 1 and
51 for both species, through task 49 handling shard numbers 49 and 99. Each case
finishes and closes its worker pool before the next starts in the same allocated
node. All 50 tasks together cover all 200 scientific cases. Each node runs
independent multiprocessing; this is not one shared multi-node Python process.

SeaWulf's [queue guide](https://rci.stonybrook.edu/docs/queue-selection-guide)
specifies at most 100 queued jobs per person. Slurm counts each array element
against submit limits; `%8` limits running tasks, not submitted tasks (see
[Slurm job arrays](https://slurm.schedmd.com/job_array.html)). Therefore the default
submission contains 50 allocation tasks, indices 0-49, which process the complete
analysis without later range submissions. Other running/pending jobs still use
submit slots; if fewer than 50 slots are available, inspect `squeue -r -u "$USER"`
and wait for enough existing jobs to end. A dependency chain of additional queued
batches would still consume submit slots and is unnecessary with this arrangement.

Protein sequence identity is the stopping criterion. The COSMIC v3.3.1 GRCh37
profile table is the same fixed table used by the checked human runner. GRCh37
identifies that mutation-profile table, not bat assembly coordinates; sequences
remain HLmolMol2 and HLmyoMyo6. All transcript isoforms are retained, including
383 Molossus / 79 Myotis records with internal stops. There is no gene/orthologue
selection, synonymous-codon null, BLOSUM criterion or experimental enzyme profile.

## 1. Update and submit the array

```bash
cd /gpfs/projects/DamaghiGroup/projects/sinabro_pub
git -C Simple_robustness_analysis_of_cds pull --ff-only
mkdir -p logs

# The pilot already passed on this project; use your own passing pilot elsewhere.
export SINABRO_BAT_PILOT="$PWD/runs/bat_cosmic_pilot_2324062"
export SINABRO_N_TRAJ=1000
export SINABRO_SEED=42

sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/run_bat_cosmic_seawulf.slurm
```

This single submission covers both species completely. **Do not submit later
50-99 / 100-149 / 150-199 ranges from the earlier instructions.** The submitted
array index now identifies an allocation task, not a single species/shard case.
run_bat_cosmic_tasks.py assigns multiple cases to each allocation while retaining
the original run_bat_cosmic.py scientific shard mapping and checkpoint logic.
The Slurm wrapper fixes task-count=50; it does not infer that value from the
number of array elements when retrying a subset. When continuing in another
shell, restore the same pilot/count/seed/run settings first. A rejected earlier
submission started no calculation and requires no output cleanup.

The driver requires a passing pilot with at least 100 paths per combination and
checks both extraction-manifest and selected-FASTA hashes. It verifies frozen
IDs, full FASTA/QC hashes, counts and lengths before calling existing run_rums.py.
Checks run on allocated compute nodes. The wrapper activates PROJECT/sinabro-revision
and uses PROJECT/sinabro; no package installation occurs.

Defaults:

```text
data/processed/bat1k_historical_v1/{molMol2,myoMyo6}/
runs/hpc_bat_cosmic_n1000_seed42/{molMol2,myoMyo6}/
  bat_input_provenance.json
  run_manifest.json
  transcripts/*.npz
  shards/000.done.json ... 099.done.json
logs/bat-cosmic-ARRAYJOB_ARRAYINDEX.{out,err}
```

Overrides: PROJECT, ANALYSIS_DIR, SINABRO_ENV, SINABRO_SOURCE,
SINABRO_BAT_PROCESSED, SINABRO_BAT_RUN, SINABRO_N_TRAJ, SINABRO_SEED,
SINABRO_MAX_STEPS and SINABRO_PYTHON. Use the same input/run/count settings for
later audit/export. Scientific changes require a documented new run.

## 2. Completion and continuation

Wait for the single 50-task array to end before auditing. Inspect logs and
scheduler state:

```bash
# Replace ARRAY_JOB_ID with the sbatch response.
squeue -j ARRAY_JOB_ID
sacct -j ARRAY_JOB_ID --format=JobID,State,ExitCode,Elapsed
```

Each species must ultimately have 100 done records and 31,185/23,213 checkpoints.
Done markers alone do not prove all paths are complete; stage 3 reads every path.
If interrupted, wait for old elements to stop, then resubmit the same command
with unchanged code, package, environment, pilot, inputs, seed, repeat/shard counts
and output root. Completed transcript files are checked/skipped; interrupted
transcripts are recomputed. Each retried allocation task visits the same four
scientific cases and skips already completed transcripts in each. Existing
checkpoints are not overwritten. Configuration
or provenance mismatches are rejected. Stored incomplete scientific paths are not
silently converted into successes by a resume; investigate before changing a cap.

For failed allocation tasks only, override indices within 0-49, e.g.
`sbatch --array=3,8%8 .../run_bat_cosmic_seawulf.slurm`. Keep both task-count=50 and
the 100-shard scientific mapping unchanged, even for a two-element retry. New scientific settings/package/driver
versions should use a new SINABRO_BAT_RUN root. Earlier outputs remain intact.

## 3. Audit every saved path

After all 50 allocation tasks complete (all 100 shards per species):

```bash
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/validate_bat_results_seawulf.slurm
```

Default output: runs/bat_cosmic_audit_JOBID/. process_bat_results.py calls the
existing read-only archive validator for each species. It checks every member,
per-path counts/terminal reasons, configuration/shard binding, file sets, CRC,
and recomputed mean/SD/median/Monte Carlo SE. No simulation runs in this stage.

Success: audit_manifest.json stage_passed=true, and both validation_summary.json
files report validation_passed=true/all_paths_complete=true. Expected complete
paths: Molossus 249,480,000; Myotis 185,704,000; total 435,184,000. Missing,
invalid, no-available-event and capped path counts must all be zero. This final
audit establishes completeness, including results from resumed invocations.

## 4. Export tables only after audit success

```bash
# Replace AUDIT_JOB_ID with the successful audit's number.
export SINABRO_BAT_AUDIT="$PWD/runs/bat_cosmic_audit_AUDIT_JOB_ID"
sbatch Simple_robustness_analysis_of_cds/06_Bat_genome/author_hpc/aggregate_bat_results_seawulf.slurm
```

Default output: runs/bat_cosmic_export_JOBID/. Per species: four canonical metric
tables, SPECIES_rums_profile.csv (historical-format means), transcript_annotations.csv,
signature_summary.csv and aggregation_manifest.json. The exporter binds every
checkpoint to its audit-time SHA256; changed/missing files or incomplete paths
prevent publication of that species' tables. Output directories must be new.
export_manifest.json gives overall two-species success; a successful species is
preserved if the other fails. SINABRO_BAT_AUDIT_OUTPUT and SINABRO_BAT_EXPORT_OUTPUT
override the new audit/export destinations respectively.

Bat mRNA IDs are retained with species, coordinates, strand, CDS length and internal
stop flags. Gene/Ensembl/orthologue IDs are not invented; gene_mapping_available=false.
Statistics equally weight transcript isoforms. Within-path SD/SE and dispersion
across transcripts are distinct. These are intermediate tables; later comparative
analysis and adaptive interpretations are not completed by export. Human values
for the same eight profiles already exist and need no new human simulation.

## Local verification before publication

Beyond the passing HPC pilot (job 2324062), all 462 internal-stop transcripts plus
the six shortest/three longest records per species were checked locally under all
eight profiles with two paths per combination: 480 records, 7,680 paths, no incomplete
paths. Both archive audits and bat-specific table exports passed. CDS lengths ranged
from 78 to 104,835 nt. This is a bounded check, not production computation.

Portable entrypoints: run_bat_cosmic_tasks.py --help (complete submission plan),
run_bat_cosmic.py --help (one scientific case), process_bat_results.py --help and
aggregate_bat_results.py --help. For bounded export checks only, --input together
with --allow-subset requires exact sequence identity to the full validated cohort
and explicitly labels a subset population. Production Slurm never selects subsets.

## The six experimental APOBEC profiles

The workbook measures 16 NNC preferences for three enzyme families from each bat:
six empirical profiles, not six COSMIC columns. The Myotis APOBEC3A-like source is
labelled APOBEC3A3 and retains that identity. Each profile can be applied to human,
Molossus and Myotis CDSs for native/foreign-enzyme and human-control comparisons.
COSMIC SBS2/SBS13 do not replace these measurements. Normalization, NNC target
positions, strand expansion and historical conversion choices need separate
validation before those runs.
