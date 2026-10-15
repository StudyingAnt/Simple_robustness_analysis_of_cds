#!/usr/bin/env bash
# Submit from the analysis repository root after activating the Python environment.
# mkdir -p logs
# sbatch 02_run_sinabro/run_rums_seawulf.slurm
# A new submission with the SAME configuration resumes completed transcripts.
#SBATCH --job-name=sinabro-rums
#SBATCH --partition=short-96core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --array=0-99%8
#SBATCH --time=04:00:00
#SBATCH --output=logs/rums-%A_%a.out
#SBATCH --error=logs/rums-%A_%a.err

set -euo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch; do not run on a login node}"
: "${SLURM_ARRAY_TASK_ID:?This script requires a Slurm array}"
: "${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
cd -- "$SLURM_SUBMIT_DIR"

# Override these with exported environment variables before sbatch if needed.
PYTHON_BIN="${SINABRO_PYTHON:-python}"
PACKAGE_SOURCE="${SINABRO_SOURCE:-$(dirname "$PWD")/sinabro}"
INPUT_FASTA="${SINABRO_FASTA:-$PWD/data/processed/gencode_v40/gencode.v40.pc_transcripts.nopary.cdsplus.fa}"
OUTPUT_ROOT="${SINABRO_OUTPUT:-$PWD/results/hpc_gencode_v40_n1000}"
N_TRAJ="${SINABRO_N_TRAJ:-1000}"
MASTER_SEED="${SINABRO_SEED:-42}"
MAX_STEPS="${SINABRO_MAX_STEPS:-10000}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID on $(hostname)"
echo "Input: $INPUT_FASTA; output: $OUTPUT_ROOT; paths per combination: $N_TRAJ"

"$PYTHON_BIN" - "$INPUT_FASTA" "$OUTPUT_ROOT" "$PACKAGE_SOURCE" "$N_TRAJ" "$MASTER_SEED" "$MAX_STEPS" "$0" <<'PY'
import csv
import fcntl
import hashlib
import heapq
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import re
import sys
import tempfile
import time

import numpy as np
import Bio
from Bio import SeqIO

source = Path(sys.argv[3]).resolve()
if not (source / 'sinabro' / '__init__.py').is_file():
    raise SystemExit(f'Set SINABRO_SOURCE to the Sinabro 3.0.0 checkout: {source}')
sys.path.insert(0, str(source))
import sinabro
from sinabro.core import signature_path, full_signature

SHARDS = 100  # Remains 100 even when resubmitting a subset of array indices.
N_TRAJ, SEED, MAX_STEPS = map(int, sys.argv[4:7])
if sinabro.__version__ != '3.0.0':
    raise SystemExit('This manuscript runner requires Sinabro 3.0.0.')
if min(N_TRAJ, MAX_STEPS) < 1 or SEED < 0:
    raise SystemExit('Invalid repeat count, seed, or safety cap.')
INPUT, OUTPUT = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
TASK = int(os.environ['SLURM_ARRAY_TASK_ID'])
WORKERS = int(os.environ['SLURM_CPUS_PER_TASK'])
if not 0 <= TASK < SHARDS or WORKERS < 1:
    raise SystemExit('Invalid Slurm array index or allocated CPUs.')
PROFILE = signature_path('3.3.1', 'GRCh37')
with PROFILE.open() as stream:
    SIGNATURES = next(csv.reader(stream, delimiter='\t'))[1:]
if len(SIGNATURES) != 79 or len(set(SIGNATURES)) != 79:
    raise SystemExit('Expected the complete 79-profile COSMIC v3.3.1 GRCh37 table.')


def digest_file(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def package_digest():
    result = hashlib.sha256()
    for path in sorted((source / 'sinabro').rglob('*.py')):
        result.update(str(path.relative_to(source)).encode() + b'\0')
        result.update(path.read_bytes())
    return result.hexdigest()


def atomic_json(path, value):
    fd, name = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


CONFIG = dict(
    runner='seawulf-rums-v1', input_sha256=digest_file(INPUT),
    runner_python_sha256=hashlib.sha256(
        Path(sys.argv[7]).read_text().split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0].encode()
    ).hexdigest(),
    signature_sha256=digest_file(PROFILE), signatures=SIGNATURES,
    n_trajectories=N_TRAJ, master_seed=SEED, max_steps=MAX_STEPS,
    cosmic_version='3.3.1', genome_ref='GRCh37',
    criterion='protein_identity', cds_start=1, cds_end='len(sequence)-1',
    translation_table=1, internal_stops_retained_as_preprocessed=True,
    sinabro=sinabro.__version__, package_source_sha256=package_digest(),
    python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__,
    shards=SHARDS, assignment='greedy_length_balance_v1',
)
FINGERPRINT = hashlib.sha256(json.dumps(CONFIG, sort_keys=True).encode()).hexdigest()


def checkpoint_path(identifier):
    return OUTPUT / 'transcripts' / (identifier + '.npz')


def check_checkpoint(path, identifier):
    if not path.exists():
        return False
    with np.load(path, allow_pickle=False) as saved:
        metadata = json.loads(str(saved['metadata'].item()))
        if metadata['config_sha256'] != FINGERPRINT or metadata['transcript_id'] != identifier:
            raise RuntimeError(f'Checkpoint configuration mismatch: {path}')
        counts, codes = saved['accepted_steps'], saved['terminal_codes']
        if counts.shape != (N_TRAJ, len(SIGNATURES)) or codes.shape != counts.shape:
            raise RuntimeError(f'Checkpoint shape mismatch: {path}')
        if not np.isin(codes, [0, 1, 2]).all():
            raise RuntimeError(f'Unfinished checkpoint: {path}')
        if ((codes == 0) & (counts < 0)).any() or ((codes != 0) & (counts != -1)).any():
            raise RuntimeError(f'Invalid checkpoint counts: {path}')
    return True


def compute_case(job):
    identifier, sequence, signature = job
    # Same case seeds and SeedSequence prefixes as benchmark_pilot.py.
    key = f'{SEED}\0{identifier}\0{signature}'
    case_seed = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], 'big')
    seeds = np.random.SeedSequence(case_seed).spawn(N_TRAJ)
    full_signature(PROFILE, signature, .5)
    counts = np.full(N_TRAJ, -1, dtype=np.int32)
    codes = np.zeros(N_TRAJ, dtype=np.uint8)
    started = time.perf_counter()
    for index, seed in enumerate(seeds):
        trajectory = sinabro.Trajectory(identifier, sequence, seed=seed,
                                       cds_start=1, cds_end=len(sequence) - 1)
        trajectory.autofill(condition='nonsynonymous', method='signature',
                            mutational_signature=signature, cosmic_version='3.3.1',
                            genome_ref='GRCh37', max_steps=MAX_STEPS)
        reason = trajectory.terminal_reason
        if reason == 'nonsynonymous':
            counts[index] = trajectory.accepted_steps
        elif reason == 'no_available_event':
            codes[index] = 1
        elif reason == 'max_steps':
            codes[index] = 2
        else:
            raise RuntimeError(f'Unexpected terminal reason: {reason}')
    return identifier, signature, counts, codes, time.perf_counter() - started


def save_transcript(identifier, fasta_id, item):
    counts, codes = item['counts'], item['codes']
    complete = (codes == 0).sum(axis=0)
    valid = complete == N_TRAJ
    mean, sd, median, se = [np.full(len(SIGNATURES), np.nan) for _ in range(4)]
    if valid.any():
        mean[valid] = counts[:, valid].mean(axis=0)
        median[valid] = np.median(counts[:, valid], axis=0)
        if N_TRAJ > 1:
            sd[valid] = counts[:, valid].std(axis=0, ddof=1)
            se[valid] = sd[valid] / np.sqrt(N_TRAJ)
    metadata = dict(config_sha256=FINGERPRINT, transcript_id=identifier, fasta_id=fasta_id,
                    n_trajectories=N_TRAJ, signatures=SIGNATURES,
                    terminal_codes={'0': 'nonsynonymous', '1': 'no_available_event', '2': 'max_steps'},
                    completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    target = checkpoint_path(identifier)
    fd, temporary = tempfile.mkstemp(prefix='.' + identifier, suffix='.tmp.npz', dir=target.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            np.savez_compressed(stream, accepted_steps=counts, terminal_codes=codes,
                                n_complete=complete, n_incomplete=N_TRAJ - complete,
                                mean=mean, sd=sd, median=median, mc_se=se,
                                simulation_seconds=item['seconds'], metadata=json.dumps(metadata))
            stream.flush()
            os.fsync(stream.fileno())
        # A hard link publishes a complete file atomically and refuses overwrite.
        os.link(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return int((codes != 0).sum())


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / 'transcripts').mkdir(exist_ok=True)
    (OUTPUT / 'shards').mkdir(exist_ok=True)
    # Serialize common-manifest creation; reject mixed code/data/settings.
    with (OUTPUT / '.configuration.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifest = OUTPUT / 'run_manifest.json'
        if manifest.exists():
            existing = json.loads(manifest.read_text())
            if existing['config_sha256'] != FINGERPRINT:
                raise RuntimeError('Output belongs to a different run. Choose a new SINABRO_OUTPUT.')
        else:
            atomic_json(manifest, dict(config_sha256=FINGERPRINT, configuration=CONFIG,
                                      input_path=str(INPUT), package_source=str(source)))
    # Prevent simultaneous duplicate submissions of this same shard.
    shard_lock = (OUTPUT / 'shards' / f'{TASK:03d}.lock').open('a+')
    fcntl.flock(shard_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    records = list(SeqIO.parse(INPUT, 'fasta'))
    identifiers = [record.id.split('|')[0] for record in records]
    if len(records) != 73214 or len(set(identifiers)) != len(records):
        raise RuntimeError('Expected 73,214 uniquely identified GENCODE v40 CDS-plus records.')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', identifier) for identifier in identifiers):
        raise RuntimeError('Unsafe transcript identifier.')
    heap = [(0, index) for index in range(SHARDS)]
    heapq.heapify(heap)
    assigned = []
    for record in sorted(records, key=lambda row: (-(len(row.seq) - 2), row.id)):
        load, shard = heapq.heappop(heap)
        if shard == TASK:
            assigned.append(record)
        heapq.heappush(heap, (load + len(record.seq) - 2, shard))
    todo = [record for record in assigned
            if not check_checkpoint(checkpoint_path(record.id.split('|')[0]), record.id.split('|')[0])]
    skipped = len(assigned) - len(todo)
    names = {record.id.split('|')[0]: record.id for record in todo}
    print(f'Shard {TASK}: assigned={len(assigned)}, resume_skip={skipped}, remaining={len(todo)}, workers={WORKERS}', flush=True)
    jobs = ((record.id.split('|')[0], str(record.seq), signature)
            for record in todo for signature in SIGNATURES)
    pending, processed, incomplete = {}, 0, 0
    signature_indices = {name: index for index, name in enumerate(SIGNATURES)}
    with mp.get_context('fork').Pool(processes=WORKERS) as pool:
        for identifier, signature, counts, codes, seconds in pool.imap_unordered(compute_case, jobs, chunksize=1):
            if identifier not in pending:
                pending[identifier] = dict(
                    counts=np.full((N_TRAJ, len(SIGNATURES)), -2, dtype=np.int32),
                    codes=np.full((N_TRAJ, len(SIGNATURES)), 255, dtype=np.uint8),
                    seconds=np.zeros(len(SIGNATURES)), received=set())
            item = pending[identifier]
            column = signature_indices[signature]
            if column in item['received']:
                raise RuntimeError('Duplicate case result.')
            item['counts'][:, column], item['codes'][:, column] = counts, codes
            item['seconds'][column] = seconds
            item['received'].add(column)
            if len(item['received']) == len(SIGNATURES):
                incomplete += save_transcript(identifier, names[identifier], item)
                del pending[identifier]
                processed += 1
                if processed % 10 == 0 or processed == len(todo):
                    print(f'Shard {TASK}: saved {processed}/{len(todo)} transcripts; incomplete paths={incomplete}', flush=True)
    if pending or processed != len(todo):
        raise RuntimeError('Shard ended with unfinished transcripts.')
    atomic_json(OUTPUT / 'shards' / f'{TASK:03d}.done.json',
                dict(config_sha256=FINGERPRINT, assigned=len(assigned), skipped=skipped,
                     newly_processed=processed, incomplete_paths_in_new_results=incomplete,
                     slurm_job_id=os.environ['SLURM_JOB_ID']))
    print(f'Shard {TASK} complete. Output: {OUTPUT}', flush=True)


if __name__ == '__main__':
    main()
PY
