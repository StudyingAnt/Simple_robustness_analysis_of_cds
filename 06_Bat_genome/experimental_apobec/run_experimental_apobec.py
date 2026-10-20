"""Simulate six experimental NNC enzyme profiles with deterministic path seeds."""
import argparse
from contextlib import ExitStack
import hashlib
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

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
import Bio
from Bio import SeqIO

from prepare_profiles import digest, load_profiles, require

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / '02_run_sinabro'))
import run_rums as storage_helpers
from run_rums import atomic_json, assign_records, exclusive_lock


def initialize_worker(settings):
    global SETTINGS, SINABRO, PROFILE_BY_ID, MOTIFS
    SETTINGS = settings
    storage_helpers.load_package(Path(settings['source']))
    SINABRO = storage_helpers.sinabro
    MOTIFS = settings['profile_body']['motifs']
    PROFILE_BY_ID = {row['id']: row['weights'] for row in settings['profile_body']['profiles']}


def compute_case(job):
    identifier, sequence, profile = job
    key = f'{SETTINGS["seed"]}\0{SETTINGS["dataset"]}\0{identifier}\0{profile}'
    seeds = np.random.SeedSequence(int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], 'big')).spawn(SETTINGS['n_traj'])
    counts = np.full(SETTINGS['n_traj'], -1, dtype=np.int32)
    codes = np.zeros(SETTINGS['n_traj'], dtype=np.uint8)
    started, epoch, cpu = time.perf_counter(), time.time(), time.process_time()
    for index, seed in enumerate(seeds):
        trajectory = SINABRO.Trajectory(identifier, sequence, seed=seed, cds_start=1, cds_end=len(sequence) - 1)
        trajectory.autofill(condition='nonsynonymous', method='mut_types',
                            mut_types=MOTIFS, mut_type_probs=PROFILE_BY_ID[profile],
                            strand='both', strand_bias=.5, max_steps=SETTINGS['max_steps'])
        if trajectory.terminal_reason == 'nonsynonymous':
            counts[index] = trajectory.accepted_steps
        elif trajectory.terminal_reason == 'no_available_event':
            codes[index] = 1
        elif trajectory.terminal_reason == 'max_steps':
            codes[index] = 2
        else:
            raise RuntimeError(f'Unexpected terminal reason: {trajectory.terminal_reason}')
    return identifier, profile, counts, codes, time.perf_counter() - started, dict(
        worker_pid=os.getpid(), started_epoch=epoch, finished_epoch=time.time(), cpu_seconds=time.process_time() - cpu)


def save_checkpoint(path, fasta_id, values, config, fingerprint):
    counts, codes = values['counts'], values['codes']
    n, size = config['n_trajectories'], len(config['signatures'])
    complete = (codes == 0).sum(axis=0)
    valid = complete == n
    stats = {name: np.full(size, np.nan) for name in ('mean', 'sd', 'median', 'mc_se')}
    if valid.any():
        selected = counts[:, valid]
        stats['mean'][valid], stats['median'][valid] = selected.mean(axis=0), np.median(selected, axis=0)
        if n > 1:
            stats['sd'][valid] = selected.std(axis=0, ddof=1)
            stats['mc_se'][valid] = stats['sd'][valid] / np.sqrt(n)
    metadata = dict(config_sha256=fingerprint, transcript_id=fasta_id.split('|')[0], fasta_id=fasta_id,
                    n_trajectories=n, signatures=config['signatures'],
                    terminal_codes={'0': 'nonsynonymous', '1': 'no_available_event', '2': 'max_steps'})
    fd, temporary = tempfile.mkstemp(prefix='.' + path.stem, suffix='.tmp.npz', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            np.savez_compressed(stream, accepted_steps=counts, terminal_codes=codes, n_complete=complete,
                                n_incomplete=n - complete, simulation_seconds=values['seconds'],
                                metadata=json.dumps(metadata), **stats)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return int((codes != 0).sum())


def check_existing(path, fasta_id, config, fingerprint):
    if not path.exists():
        return False
    # Reuse the independent archive validator instead of accepting filenames alone.
    from validate_results import audit_checkpoint
    result = audit_checkpoint((path, fasta_id), config, fingerprint)
    require(result['structurally_valid'], f'Invalid existing checkpoint: {path}: {result["error"]}')
    return True


def run(args):
    require(min(args.workers, args.n_traj, args.max_steps, args.n_shards) > 0 and args.seed >= 0,
            'Invalid worker/path/shard/cap/seed settings')
    require(0 <= args.shard_index < args.n_shards, 'Shard index outside range')
    body = load_profiles(args.profiles_dir)
    profiles = [row['id'] for row in body['profiles']]
    storage_helpers.load_package(args.sinabro_source)
    source = storage_helpers.source
    input_path, output = args.input.resolve(), args.out_dir.resolve()
    require(output != args.profiles_dir.resolve() and args.profiles_dir.resolve() not in output.parents,
            'Run output must be outside prepared profiles')
    require(output != input_path.parent and input_path.parent not in output.parents,
            'Run output must be outside the input directory')
    config = dict(runner='rums-v3', mutation_model='experimental_apobec_NNC_v1', dataset=args.dataset,
                  input_sha256=digest(input_path), runner_python_sha256=digest(__file__),
                  profiles_sha256=digest(args.profiles_dir / 'profiles.json'),
                  profile_preparation_sha256=digest(args.profiles_dir / 'preparation_manifest.json'),
                  profile_loader_sha256=digest(HERE / 'prepare_profiles.py'),
                  signatures=profiles, n_trajectories=args.n_traj, master_seed=args.seed,
                  seed_key='master_seed,dataset,transcript_id,enzyme_id; sha256-first8-big; SeedSequence prefixes',
                  max_steps=args.max_steps, criterion='protein_identity', cds_start=1, cds_end='len(sequence)-1',
                  translation_table=1, internal_stops_retained_as_preprocessed=True,
                  motif_weighting=body['weighting'], context_policy=body['context_policy'],
                  strand='both', strand_bias=.5, cosmic_version=None, genome_ref=None,
                  python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__,
                  sinabro=storage_helpers.sinabro.__version__, package_source_sha256=storage_helpers.package_digest(),
                  shards=args.n_shards, assignment='greedy_length_balance_v1')
    fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    records = list(SeqIO.parse(input_path, 'fasta'))
    ids = [r.id.split('|')[0] for r in records]
    require(records and len(ids) == len(set(ids)) and all(re.fullmatch('[A-Za-z0-9_.-]+', i) for i in ids),
            'Empty/duplicate/unsafe transcript IDs')
    require(args.expected_records is None or len(records) == args.expected_records, 'Unexpected input count')
    output.mkdir(parents=True, exist_ok=True)
    for name in ('transcripts', 'shards'):
        (output / name).mkdir(exist_ok=True)
    with exclusive_lock(output / '.configuration.lock'):
        manifest = output / 'run_manifest.json'
        if manifest.exists():
            require(json.loads(manifest.read_text())['config_sha256'] == fingerprint,
                    'Output configuration changed; choose a new run directory')
        else:
            atomic_json(manifest, dict(config_sha256=fingerprint, configuration=config,
                                      input_path=str(input_path), package_source=str(source)))
    settings = dict(source=str(source), dataset=args.dataset, n_traj=args.n_traj,
                    seed=args.seed, max_steps=args.max_steps, profile_body=body)
    initialize_worker(settings)
    started = time.perf_counter()
    with exclusive_lock(output / 'shards' / f'{args.shard_index:03d}.lock', blocking=False):
        assigned = assign_records(records, args.n_shards, args.shard_index)
        todo = [r for r in assigned if not check_existing(output / 'transcripts' / (r.id.split('|')[0] + '.npz'),
                                                         r.id, config, fingerprint)]
        skipped = len(assigned) - len(todo)
        print(f'{args.dataset} shard {args.shard_index}: assigned={len(assigned)}, '
              f'resume_skip={skipped}, remaining={len(todo)}, workers={args.workers}', flush=True)
        names = {r.id.split('|')[0]: r.id for r in todo}
        jobs = ((r.id.split('|')[0], str(r.seq), profile) for r in todo for profile in profiles)
        pending, processed, incomplete = {}, 0, 0
        metrics_path = output / 'shards' / f'{args.shard_index:03d}.cases.{uuid.uuid4().hex}.jsonl'
        simulation_start = time.perf_counter()
        with ExitStack() as stack:
            metrics = stack.enter_context(metrics_path.open('x'))
            pool = stack.enter_context(mp.get_context('spawn').Pool(args.workers, initializer=initialize_worker,
                                                                  initargs=(settings,))) if args.workers > 1 and todo else None
            results = pool.imap_unordered(compute_case, jobs, chunksize=1) if pool else map(compute_case, jobs)
            for identifier, profile, counts, codes, seconds, execution in results:
                metrics.write(json.dumps(dict(transcript_id=identifier, signature=profile,
                                              simulation_seconds=seconds, **execution)) + '\n')
                metrics.flush()
                if identifier not in pending:
                    pending[identifier] = dict(counts=np.full((args.n_traj, 6), -2, dtype=np.int32),
                        codes=np.full((args.n_traj, 6), 255, dtype=np.uint8), seconds=np.zeros(6), received=set())
                item = pending[identifier]
                column = profiles.index(profile)
                require(column not in item['received'], 'Duplicate simulation case')
                item['counts'][:, column], item['codes'][:, column] = counts, codes
                item['seconds'][column] = seconds
                item['received'].add(column)
                if len(item['received']) == 6:
                    incomplete += save_checkpoint(output / 'transcripts' / (identifier + '.npz'), names[identifier],
                                                  item, config, fingerprint)
                    del pending[identifier]
                    processed += 1
                    if processed % 10 == 0 or processed == len(todo):
                        print(f'{args.dataset} shard {args.shard_index}: saved {processed}/{len(todo)}; '
                              f'incomplete paths={incomplete}', flush=True)
        require(not pending and processed == len(todo), 'Unfinished transcripts at shard exit')
        require(digest(input_path) == config['input_sha256']
                and digest(args.profiles_dir / 'profiles.json') == config['profiles_sha256'],
                'Input/profile changed during simulation')
        atomic_json(output / 'shards' / f'{args.shard_index:03d}.done.json', dict(
            config_sha256=fingerprint, assigned=len(assigned), skipped=skipped, newly_processed=processed,
            incomplete_paths_in_new_results=incomplete, workers=args.workers, start_method='spawn',
            elapsed_seconds=time.perf_counter() - started,
            compute_and_save_seconds=time.perf_counter() - simulation_start,
            case_metrics_file=metrics_path.name, host=platform.node()))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--profiles-dir', type=Path, required=True)
    parser.add_argument('--dataset', choices=('human', 'molMol2', 'myoMyo6'), required=True)
    parser.add_argument('--sinabro-source', type=Path)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-shards', type=int, default=1)
    parser.add_argument('--shard-index', type=int, default=0)
    parser.add_argument('--expected-records', type=int)
    run(parser.parse_args(argv))


if __name__ == '__main__':
    mp.freeze_support()
    main()
