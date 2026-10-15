"""Reproduce RUMS on a local computer or a scheduler-managed compute node."""
import argparse
from contextlib import contextmanager, ExitStack
import csv
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
import uuid

# Each process performs serial numerical work; multiprocessing controls concurrency.
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
import Bio
from Bio import SeqIO


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


@contextmanager
def exclusive_lock(path, blocking=True):
    """Lock an open file on POSIX or Windows; the OS releases locks after a crash."""
    with Path(path).open('a+b') as handle:
        if os.name == 'nt':
            import msvcrt
            if handle.seek(0, os.SEEK_END) == 0:
                handle.write(b'\0')
                handle.flush()
            while True:
                handle.seek(0)
                try:
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError:
                    if not blocking:
                        raise
                    time.sleep(0.1)
            try:
                yield
            finally:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            operation = fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB)
            fcntl.flock(handle, operation)
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Preprocessed CDS-plus FASTA')
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--sinabro-source', type=Path, help='Optional package checkout; otherwise use installed Sinabro')
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-shards', type=int, default=1)
    parser.add_argument('--shard-index', type=int, default=0)
    parser.add_argument('--expected-records', type=int, help='Optional input-record count assertion')
    parser.add_argument('--signatures', nargs='+', help='Default: all 79 COSMIC v3.3.1 GRCh37 profiles')
    parser.add_argument('--start-method', choices=mp.get_all_start_methods(), default='spawn')
    args = parser.parse_args(argv)
    if min(args.n_traj, args.max_steps, args.workers, args.n_shards) < 1 or args.seed < 0:
        parser.error('Repeat count, safety cap, workers and shard count must be positive; seed must be nonnegative.')
    if not 0 <= args.shard_index < args.n_shards:
        parser.error('Shard index must be between 0 and n-shards minus one.')
    if args.expected_records is not None and args.expected_records < 1:
        parser.error('Expected record count must be positive.')
    if args.signatures is not None and len(set(args.signatures)) != len(args.signatures):
        parser.error('Signature names must be unique.')
    return args


def load_package(checkout):
    global source, sinabro, full_signature, signature_path
    if checkout is not None:
        checkout = Path(checkout).resolve()
        if not (checkout / 'sinabro' / '__init__.py').is_file():
            raise ValueError(f'Invalid --sinabro-source checkout: {checkout}')
        sys.path.insert(0, str(checkout))
    import sinabro
    from sinabro.core import full_signature, signature_path
    if getattr(sinabro, '__version__', None) != '3.0.0' or sinabro.__file__ is None:
        raise ValueError('This manuscript runner requires an importable Sinabro 3.0.0 package.')
    source = Path(sinabro.__file__).resolve().parent.parent


def initialize_worker(settings):
    """Initialize spawned workers without reading/hashing the full FASTA again."""
    global N_TRAJ, SEED, MAX_STEPS, PROFILE
    load_package(settings['source'])
    N_TRAJ, SEED, MAX_STEPS = settings['n_traj'], settings['seed'], settings['max_steps']
    PROFILE = Path(settings['profile'])


def configure(args):
    global N_TRAJ, SEED, MAX_STEPS, INPUT, OUTPUT, TASK, WORKERS, SHARDS, START_METHOD
    global PROFILE, SIGNATURES, CONFIG, FINGERPRINT
    load_package(args.sinabro_source)
    N_TRAJ, SEED, MAX_STEPS = args.n_traj, args.seed, args.max_steps
    INPUT, OUTPUT = args.input.resolve(), args.out_dir.resolve()
    TASK, WORKERS, SHARDS = args.shard_index, args.workers, args.n_shards
    START_METHOD = args.start_method
    PROFILE = signature_path('3.3.1', 'GRCh37')
    with PROFILE.open() as stream:
        available = next(csv.reader(stream, delimiter='\t'))[1:]
    if len(available) != 79 or len(set(available)) != 79:
        raise ValueError('Expected the complete 79-profile COSMIC v3.3.1 GRCh37 table.')
    SIGNATURES = args.signatures or available
    if set(SIGNATURES) - set(available):
        raise ValueError(f'Unknown signatures: {sorted(set(SIGNATURES) - set(available))}')

    CONFIG = dict(
        runner='rums-v3', input_sha256=digest_file(INPUT),
        runner_python_sha256=digest_file(__file__),
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


def assign_records(records, n_shards, shard_index):
    """Greedy length balancing with the original deterministic tie breaking."""
    heap = [(0, index) for index in range(n_shards)]
    heapq.heapify(heap)
    assigned = []
    for record in sorted(records, key=lambda row: (-(len(row.seq) - 2), row.id)):
        load, shard = heapq.heappop(heap)
        if shard == shard_index:
            assigned.append(record)
        heapq.heappush(heap, (load + len(record.seq) - 2, shard))
    return assigned


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
    started_epoch, cpu_started = time.time(), time.process_time()
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
    return (identifier, signature, counts, codes, time.perf_counter() - started,
            dict(worker_pid=os.getpid(), started_epoch=started_epoch, finished_epoch=time.time(),
                 cpu_seconds=time.process_time() - cpu_started))


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


def run_shard(args, started):
    records = list(SeqIO.parse(INPUT, 'fasta'))
    identifiers = [record.id.split('|')[0] for record in records]
    if not records or len(set(identifiers)) != len(records):
        raise ValueError('Input must contain at least one uniquely identified CDS-plus record.')
    if args.expected_records is not None and len(records) != args.expected_records:
        raise ValueError(f'Expected {args.expected_records} records; found {len(records)}.')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', identifier) for identifier in identifiers):
        raise ValueError('Unsafe transcript identifier.')
    assigned = assign_records(records, SHARDS, TASK)
    todo = [record for record in assigned
            if not check_checkpoint(checkpoint_path(record.id.split('|')[0]), record.id.split('|')[0])]
    skipped = len(assigned) - len(todo)
    names = {record.id.split('|')[0]: record.id for record in todo}
    print(f'Shard {TASK}: assigned={len(assigned)}, resume_skip={skipped}, remaining={len(todo)}, workers={WORKERS}', flush=True)
    jobs = ((record.id.split('|')[0], str(record.seq), signature)
            for record in todo for signature in SIGNATURES)
    pending, processed, incomplete = {}, 0, 0
    signature_indices = {name: index for index, name in enumerate(SIGNATURES)}
    metrics_path = OUTPUT / 'shards' / f'{TASK:03d}.cases.{uuid.uuid4().hex}.jsonl'
    simulation_start = time.perf_counter()
    with ExitStack() as stack:
        metrics = stack.enter_context(metrics_path.open('x'))
        if WORKERS == 1:
            results = map(compute_case, jobs)
        else:
            settings = dict(source=str(source), n_traj=N_TRAJ, seed=SEED,
                            max_steps=MAX_STEPS, profile=str(PROFILE))
            pool = stack.enter_context(mp.get_context(START_METHOD).Pool(
                processes=WORKERS, initializer=initialize_worker, initargs=(settings,)))
            results = pool.imap_unordered(compute_case, jobs, chunksize=1)
        for identifier, signature, counts, codes, seconds, execution in results:
            metrics.write(json.dumps(dict(transcript_id=identifier, signature=signature,
                                          simulation_seconds=seconds, **execution)) + '\n')
            metrics.flush()
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
                     workers=WORKERS, start_method=START_METHOD, host=platform.node(),
                     elapsed_seconds=time.perf_counter() - started,
                     compute_and_save_seconds=time.perf_counter() - simulation_start,
                     case_metrics_file=metrics_path.name))
    print(f'Shard {TASK} complete. Output: {OUTPUT}', flush=True)


def main(argv=None):
    args = parse_arguments(argv)
    started = time.perf_counter()
    configure(args)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / 'transcripts').mkdir(exist_ok=True)
    (OUTPUT / 'shards').mkdir(exist_ok=True)
    with exclusive_lock(OUTPUT / '.configuration.lock'):
        manifest = OUTPUT / 'run_manifest.json'
        if manifest.exists():
            existing = json.loads(manifest.read_text())
            if existing['config_sha256'] != FINGERPRINT:
                raise RuntimeError('Output belongs to a different run. Choose a new --out-dir.')
        else:
            atomic_json(manifest, dict(config_sha256=FINGERPRINT, configuration=CONFIG,
                                      input_path=str(INPUT), package_source=str(source)))
    with exclusive_lock(OUTPUT / 'shards' / f'{TASK:03d}.lock', blocking=False):
        run_shard(args, started)


if __name__ == '__main__':
    mp.freeze_support()
    main()
