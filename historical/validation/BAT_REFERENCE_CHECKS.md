# Local validation on actual Bat1K references

Date: 2026-10-17. This is a bounded execution check on a local WSL computer,
not a completed whole-genome bat simulation or a measurement on SeaWulf.
Scientific environment: Python 3.9.15, NumPy 1.21.6, Biopython 1.79, Sinabro 3.0.0.
No package or existing human-runner implementation was modified.

## Raw references and full CDS extraction

Both complete official compressed genomes, both allAnnotatedIsoforms BEDs and
the source README were downloaded and checked. Genome SHA256s are recorded in
bat1k_pilot_sources.json and pinned for subsequent downloads. They are measured
from that download, not publisher-supplied checksums. Full gzip integrity passed.

| Species | Frozen historical IDs | Extracted and valid | Required-check failures | Internal-stop transcripts retained |
| --- | ---: | ---: | ---: | ---: |
| Molossus molossus | 31,185 | 31,185 | 0 | 383 |
| Myotis myotis | 23,213 | 23,213 | 0 | 79 |

Every frozen ID has equal transcript/CDS outer coordinates, complete genomic flanks,
A/C/G/T sequence, a valid reading frame, ATG start and terminal stop. Internal stops
are retained and flagged to preserve the historical policy. No extra representative
isoform or orthologue selection was introduced.

An independent literal reconstruction of the old all-exon extraction/filter logic
was also applied to the entire annotation. Molossus accepted exactly the historical
31,185 IDs. Myotis accepted 23,215, including two additional N-containing CDSs,
mRNA9113 and mRNA9159. They are absent from the historical result table and remain
outside this reproduction cohort; why the historical results omitted them is not
recorded. All 54,398 frozen-cohort sequences agree with the reconstructed legacy
extraction after uppercasing the flanks. The Molossus original script left mixed-case
flanks in 2,638 accepted records. Exact archived historical FASTAs are unavailable,
so this does not establish their byte identity.

## Bounded COSMIC pilot

For each species, selected 12 evenly spaced CDS-length quantile representatives
without internal stops. Used SBS2, SBS4 and SBS13, seed 42, 20 trajectories per
combination and max_steps=10,000. Compared one worker with four spawned workers.
That is 720 trajectories per species per mode, 2,880 including both species/modes.

| Check | Molossus | Myotis |
| --- | --- | --- |
| Selected CDS lengths | 363-5,028 nt | 366-4,782 nt |
| Per-path accepted counts and terminal codes | Serial/parallel identical | Serial/parallel identical |
| Completed paths | All | All |
| Recomputed statistics and checkpoint archive audit | Passed in both modes | Passed in both modes |
| Observed peak simultaneous cases | 4 | 4 |

The first measured invocation times were 1.62 / 0.77 seconds for Molossus
(serial / four workers) and 1.77 / 0.87 seconds for Myotis. These small local
measurements do not estimate full-genome runtime, cluster runtime or scientific
agreement with the old random trajectories. Full processed FASTAs retain the
internal-stop records omitted from this execution demonstration.

The HPC pilot defaults to 100 trajectories per combination and eight workers;
it has not been submitted or tested on the cluster by the assistant. Follow
HPC_REPRODUCTION.md one stage at a time. Full eight-signature and experimental
APOBEC profile runs are subsequent stages.

## Automated checks

All 71 analysis tests passed, including 13 new checks for spliced plus/minus
coordinates and flanks, cohort/source integrity, ambiguous sequences, internal
stops, missing contigs/flanks, corrupt gzip data, immutable output directories,
length-based pilot selection and changed per-path result detection. All three
bat Slurm wrappers passed bash -n. Existing original scripts remain unchanged.

Large reference FASTAs, processed CDSs and generated pilot archives remain outside
Git. Versioned code, cohort ID lists and instructions are sufficient to regenerate
them on another computer with the same references and environment.

## Production preparation checks, 2026-10-19

The full eight-profile launch/audit/export workflow is documented in
[COSMIC_PRODUCTION.md](../../06_Bat_genome/COSMIC_PRODUCTION.md). New automated checks bring the
analysis suite to 84 passing tests; they cover species/shard array mapping,
pilot/input provenance, scientific configuration guards, resume mismatches,
bat-specific audited export and changed-file rejection. All three new Slurm
wrappers passed bash -n. The human runner and exporter remain unchanged.

All 462 internal-stop transcripts plus short/long boundary examples were checked
with two paths under each of the eight COSMIC profiles: 392 Molossus and 88 Myotis
records, 7,680 complete paths, both archive audits and species-specific table
exports passed. These bounded exports explicitly label their subset population.
A local 100-path serial/four-worker pilot passed for the final extraction inputs.
The production driver also dispatched one real shard per species with reduced
path counts; resubmitting the completed Molossus shard skipped its checkpoint.
No full bat production run or cluster submission was performed by the assistant.

The complete submission now uses run_bat_cosmic_tasks.py: 50 allocation tasks
each process four scientific species/shard cases sequentially. All 200 original
cases are covered once, retaining the existing per-species 100-shard assignment,
seeds, provenance and checkpoint paths. Four additional allocation-plan checks
bring the suite to 88 passing tests. A real-data four-case local orchestration
check skipped two previously completed cases and completed two more small cases
with no incomplete paths. The Slurm wrapper passed Bash syntax checking.
Manual later-range submissions are no longer part of this workflow.
