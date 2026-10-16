"""Audit saved RUMS results without importing Sinabro or running simulations."""
import argparse
import csv
import functools
import hashlib
import heapq
import io
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
import Bio
from Bio import SeqIO


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def audit_checkpoint(task, configuration, fingerprint):
    """Read every archive member (including ZIP CRC checks) and recompute summaries."""
    path, fasta_id = task
    identifier = fasta_id.split('|')[0]
    result = dict(transcript_id=identifier, sha256='', structurally_valid=False,
                  complete_paths=0, no_available_event_paths=0, max_steps_paths=0, error='')
    try:
        raw = path.read_bytes()
        result['sha256'] = hashlib.sha256(raw).hexdigest()
        with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
            names = {'accepted_steps', 'terminal_codes', 'n_complete', 'n_incomplete',
                     'mean', 'sd', 'median', 'mc_se', 'simulation_seconds', 'metadata'}
            require(len(archive.files) == len(names) and set(archive.files) == names,
                    'Unexpected or missing archive members')
            saved = {name: archive[name] for name in archive.files}
        metadata = json.loads(str(saved['metadata'].item()))
        require(metadata['config_sha256'] == fingerprint, 'Configuration fingerprint mismatch')
        require(metadata['transcript_id'] == identifier, 'Transcript identifier mismatch')
        require(metadata['fasta_id'] == fasta_id, 'FASTA identifier mismatch')
        n = configuration['n_trajectories']
        signatures = configuration['signatures']
        require(metadata['n_trajectories'] == n, 'Trajectory metadata mismatch')
        require(metadata['signatures'] == signatures, 'Signature order mismatch')
        require(metadata['terminal_codes'] == {'0': 'nonsynonymous',
                    '1': 'no_available_event', '2': 'max_steps'}, 'Terminal-code metadata mismatch')
        counts, codes = saved['accepted_steps'], saved['terminal_codes']
        require(counts.shape == codes.shape == (n, len(signatures)), 'Trajectory array shape mismatch')
        require(counts.dtype == np.dtype('int32'), 'accepted_steps must be int32')
        require(codes.dtype == np.dtype('uint8'), 'terminal_codes must be uint8')
        require(np.isin(codes, [0, 1, 2]).all(), 'Unknown or unfinished terminal code')
        require(not ((codes == 0) & ((counts < 0) | (counts >= configuration['max_steps']))).any(),
                'Completed path has an invalid accepted-step count')
        require(not ((codes != 0) & (counts != -1)).any(), 'Incomplete path must have count -1')
        complete = (codes == 0).sum(axis=0)
        for name, expected in [('n_complete', complete), ('n_incomplete', n - complete)]:
            require(saved[name].dtype.kind in 'iu' and saved[name].shape == (len(signatures),)
                    and np.array_equal(saved[name], expected), name + ' mismatch')
        valid = complete == n
        expected = {name: np.full(len(signatures), np.nan)
                    for name in ('mean', 'sd', 'median', 'mc_se')}
        if valid.any():
            values = counts[:, valid]
            expected['mean'][valid] = values.mean(axis=0)
            expected['median'][valid] = np.median(values, axis=0)
            if n > 1:
                expected['sd'][valid] = values.std(axis=0, ddof=1)
                expected['mc_se'][valid] = expected['sd'][valid] / np.sqrt(n)
        for name, values in expected.items():
            actual = saved[name]
            require(actual.dtype.kind == 'f' and actual.shape == values.shape,
                    name + ' shape or dtype mismatch')
            require(np.allclose(actual, values, rtol=1e-10, atol=1e-12, equal_nan=True),
                    name + ' does not match recomputed statistics')
        seconds = saved['simulation_seconds']
        require(seconds.shape == (len(signatures),) and seconds.dtype.kind == 'f'
                and np.isfinite(seconds).all() and (seconds >= 0).all(), 'Invalid simulation timings')
        result.update(structurally_valid=True, complete_paths=int(complete.sum()),
                      no_available_event_paths=int((codes == 1).sum()),
                      max_steps_paths=int((codes == 2).sum()))
    except Exception as error:
        result['error'] = f'{type(error).__name__}: {error}'
    return result


def validate_run(run_dir, input_path, report_dir, workers=1, expected_records=None):
    """Produce reports; never write into the run directory or change checkpoints."""
    started = time.perf_counter()
    run_dir, report_dir = Path(run_dir).resolve(), Path(report_dir).resolve()
    require(report_dir != run_dir and run_dir not in report_dir.parents,
            'Report directory must be outside the simulation run directory')
    # A new directory prevents overwriting an earlier audit, even after failure.
    report_dir.mkdir(parents=True, exist_ok=False)
    summary = dict(validation_passed=False, run_directory=str(run_dir),
                   validator_sha256=digest_file(__file__), workers=workers,
                   validation_python=platform.python_version(), validation_numpy=np.__version__,
                   validation_biopython=Bio.__version__,
                   started_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                   issues=[])
    issues = summary['issues']
    try:
        require(workers > 0, 'Worker count must be positive')
        manifest_path = run_dir / 'run_manifest.json'
        manifest_bytes = manifest_path.read_bytes()
        manifest = json.loads(manifest_bytes)
        config, fingerprint = manifest['configuration'], manifest['config_sha256']
        require(hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest() == fingerprint,
                'Manifest configuration fingerprint is invalid')
        require(config['runner'] == 'rums-v3', 'Unsupported runner format')
        n, signatures, n_shards = config['n_trajectories'], config['signatures'], config['shards']
        require(n > 0 and n_shards > 0 and len(signatures) > 0
                and len(set(signatures)) == len(signatures), 'Invalid configuration dimensions')
        require(config['assignment'] == 'greedy_length_balance_v1', 'Unsupported shard assignment')
        input_path = Path(input_path or manifest['input_path']).resolve()
        input_sha = digest_file(input_path)
        require(input_sha == config['input_sha256'], 'Input FASTA SHA256 differs from the run')
        records = [(record.id, len(record.seq) - 2) for record in SeqIO.parse(input_path, 'fasta')]
        identifiers = [name.split('|')[0] for name, _ in records]
        require(records and len(set(identifiers)) == len(records), 'Empty input or duplicate transcript IDs')
        if expected_records is not None:
            require(len(records) == expected_records, 'Unexpected input transcript count')
        expected_files = {identifier + '.npz' for identifier in identifiers}
        actual_files = {path.name for path in (run_dir / 'transcripts').glob('*.npz')}
        missing, extra = sorted(expected_files - actual_files), sorted(actual_files - expected_files)
        if missing or extra:
            issues.append('Transcript file set differs from the input FASTA; see file_set.json')
        (report_dir / 'file_set.json').write_text(json.dumps(dict(missing=missing, extra=extra), indent=2) + '\n')
        heap = [(0, index) for index in range(n_shards)]
        heapq.heapify(heap)
        assigned, shard_of = [0] * n_shards, {}
        for name, length in sorted(records, key=lambda row: (-row[1], row[0])):
            load, index = heapq.heappop(heap)
            assigned[index] += 1
            shard_of[name.split('|')[0]] = index
            heapq.heappush(heap, (load + length, index))
        done = {}
        for index in range(n_shards):
            try:
                record = json.loads((run_dir / 'shards' / f'{index:03d}.done.json').read_text())
                require(record['config_sha256'] == fingerprint, 'Configuration fingerprint mismatch')
                require(record['assigned'] == assigned[index], 'Assigned transcript count mismatch')
                require(record['skipped'] >= 0 and record['newly_processed'] >= 0
                        and record['skipped'] + record['newly_processed'] == assigned[index],
                        'Processed/skipped counts mismatch')
                require(record['incomplete_paths_in_new_results'] >= 0, 'Negative incomplete count')
                done[index] = record
            except Exception as error:
                issues.append(f'Shard {index}: {type(error).__name__}: {error}')
        expected_done = {f'{i:03d}.done.json' for i in range(n_shards)}
        if {p.name for p in (run_dir / 'shards').glob('*.done.json')} - expected_done:
            issues.append('Unexpected shard completion files')
        tasks = ((run_dir / 'transcripts' / (name.split('|')[0] + '.npz'), name)
                 for name, _ in records)
        audit = functools.partial(audit_checkpoint, configuration=config, fingerprint=fingerprint)
        valid_files = invalid_files = complete = no_event = capped = 0
        incomplete_by_shard = [0] * n_shards
        fields = ['transcript_id', 'sha256', 'structurally_valid', 'complete_paths',
                  'no_available_event_paths', 'max_steps_paths', 'error']
        pool = mp.get_context('spawn').Pool(workers) if workers > 1 else None
        try:
            results = pool.imap(audit, tasks, chunksize=8) if pool else map(audit, tasks)
            with (report_dir / 'transcript_audit.csv').open('x', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                for checked, result in enumerate(results, 1):
                    writer.writerow(result)
                    valid_files += int(result['structurally_valid'])
                    invalid_files += int(not result['structurally_valid'])
                    complete += result['complete_paths']
                    no_event += result['no_available_event_paths']
                    capped += result['max_steps_paths']
                    incomplete_by_shard[shard_of[result['transcript_id']]] += (
                        result['no_available_event_paths'] + result['max_steps_paths'])
                    if checked % 1000 == 0 or checked == len(records):
                        stream.flush()
                        print(f'Checked {checked}/{len(records)}; invalid files={invalid_files}; '
                              f'incomplete paths={no_event + capped}', flush=True)
        finally:
            if pool:
                pool.close()
                pool.join()
        if invalid_files:
            issues.append('Unreadable or inconsistent transcript files; see transcript_audit.csv')
        if no_event + capped:
            issues.append('Incomplete trajectories exist; inspect terminal-code totals')
        if invalid_files == 0:
            for index, record in done.items():
                if record['skipped'] == 0 and record['incomplete_paths_in_new_results'] != incomplete_by_shard[index]:
                    issues.append(f'Shard {index}: incomplete count differs from saved arrays')
        require(manifest_path.read_bytes() == manifest_bytes, 'Manifest changed during validation')
        summary.update(config_sha256=fingerprint, configuration=config,
                       manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                       input_path=str(input_path), input_sha256=input_sha,
                       expected_transcripts=len(records), observed_transcript_files=len(actual_files),
                       missing_transcripts=len(missing), extra_transcripts=len(extra),
                       completed_shards=len(done), expected_shards=n_shards,
                       structurally_valid_files=valid_files, invalid_files=invalid_files,
                       expected_paths=len(records) * n * len(signatures),
                       checked_paths=complete + no_event + capped,
                       complete_paths=complete, no_available_event_paths=no_event,
                       max_steps_paths=capped, all_paths_complete=(invalid_files == 0 and no_event + capped == 0))
        summary['validation_passed'] = not issues
    except Exception as error:
        issues.append(f'{type(error).__name__}: {error}')
    summary['elapsed_seconds'] = time.perf_counter() - started
    summary['finished_at_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    with (report_dir / 'validation_summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--report-dir', type=Path, required=True, help='New directory outside --run-dir')
    parser.add_argument('--input', type=Path, help='Default: FASTA path recorded in run_manifest.json')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--expected-records', type=int)
    args = parser.parse_args(argv)
    if args.workers < 1 or (args.expected_records is not None and args.expected_records < 1):
        parser.error('Worker and expected-record counts must be positive')
    summary = validate_run(args.run_dir, args.input, args.report_dir, args.workers, args.expected_records)
    print(json.dumps(summary, indent=2))
    print('Report:', args.report_dir / 'validation_summary.json')
    return 0 if summary['validation_passed'] else 1


if __name__ == '__main__':
    mp.freeze_support()
    sys.exit(main())
