"""Compare one and multiple workers using the same portable RUMS runner."""
import argparse
import csv
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_rums import digest_file
from pilot_selection import SIGNATURES, select_transcripts
import numpy as np
from Bio import SeqIO


def summarize_execution(directory):
    done = json.loads((directory / 'shards/000.done.json').read_text())
    rows = [json.loads(line) for line in
            (directory / 'shards' / done['case_metrics_file']).read_text().splitlines()]
    events = []
    for row in rows:
        events.extend([(row['started_epoch'], 1), (row['finished_epoch'], -1)])
    active, peak = 0, 0
    for _, delta in sorted(events):
        active += delta
        peak = max(peak, active)
    return dict(elapsed_seconds=done['elapsed_seconds'],
                compute_and_save_seconds=done['compute_and_save_seconds'],
                workers=done['workers'], start_method=done['start_method'],
                worker_pids=sorted({row['worker_pid'] for row in rows}),
                peak_simultaneous_cases=peak,
                incomplete_paths=done['incomplete_paths_in_new_results'])


def compare_results(serial, parallel, expected_records, n_traj):
    left = {p.name: p for p in (serial / 'transcripts').glob('*.npz')}
    right = {p.name: p for p in (parallel / 'transcripts').glob('*.npz')}
    if len(left) != expected_records or set(left) != set(right):
        raise AssertionError('Missing or different transcript checkpoints.')
    checked = 0
    for name in sorted(left):
        with np.load(left[name], allow_pickle=False) as a, np.load(right[name], allow_pickle=False) as b:
            for key in ('accepted_steps', 'terminal_codes'):
                if a[key].shape != (n_traj, len(SIGNATURES)):
                    raise AssertionError(f'Unexpected shape: {name} {key}')
                np.testing.assert_array_equal(a[key], b[key], err_msg=f'{name} {key}')
            ma, mb = [json.loads(str(item['metadata'].item())) for item in (a, b)]
            if ma['signatures'] != SIGNATURES or mb['signatures'] != SIGNATURES:
                raise AssertionError('Signature order differs.')
            if ma['config_sha256'] != mb['config_sha256']:
                raise AssertionError('Scientific configurations differ.')
            checked += a['accepted_steps'].size
    return checked


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Full preprocessed GENCODE v40 CDS-plus FASTA')
    parser.add_argument('--out-dir', type=Path, required=True, help='New directory; must not exist')
    parser.add_argument('--sinabro-source', type=Path)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--n-traj', type=int, default=100)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--timeout', type=int, default=1800, help='Seconds per serial/parallel invocation')
    args = parser.parse_args(argv)
    if args.workers < 2 or min(args.n_traj, args.max_steps, args.timeout) < 1 or args.seed < 0:
        parser.error('Use at least two workers, positive repeat/cap/timeout values and a nonnegative seed.')

    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    selected, records, populations, exclusions = select_transcripts(args.input)
    selected_file = output / 'selected_transcripts.fa'
    SeqIO.write(records, selected_file, 'fasta')
    with (output / 'selected_transcripts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    runner = Path(__file__).resolve().parents[1] / 'run_rums.py'
    manifest = dict(input_path=str(args.input.resolve()), input_sha256=digest_file(args.input),
                    selected_fasta_sha256=digest_file(selected_file), selection=selected,
                    population_by_bin=populations, internal_stop_excluded_from_pilot_by_bin=exclusions,
                    signatures=SIGNATURES, n_traj_per_case=args.n_traj, master_seed=args.seed,
                    max_steps=args.max_steps, parallel_workers=args.workers, start_method='spawn',
                    validation_script_sha256=digest_file(__file__),
                    selection_script_sha256=digest_file(Path(__file__).with_name('pilot_selection.py')),
                    runner_sha256=digest_file(runner), platform=platform.platform(),
                    planned_trajectories_per_mode=len(records) * len(SIGNATURES) * args.n_traj)
    (output / 'validation_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    execution = {}
    for label, workers in [('serial', 1), (f'parallel_{args.workers}', args.workers)]:
        directory = output / label
        command = [sys.executable, str(runner), '--input', str(selected_file),
                   '--out-dir', str(directory), '--n-traj', str(args.n_traj),
                   '--seed', str(args.seed), '--max-steps', str(args.max_steps),
                   '--workers', str(workers), '--n-shards', '1', '--shard-index', '0',
                   '--expected-records', '12', '--start-method', 'spawn',
                   '--signatures', *SIGNATURES]
        if args.sinabro_source is not None:
            command.extend(['--sinabro-source', str(args.sinabro_source.resolve())])
        print(f'Running {label}: 12 transcripts x 3 signatures x {args.n_traj} paths', flush=True)
        started = time.perf_counter()
        subprocess.run(command, check=True, timeout=args.timeout)
        execution[label] = dict(summarize_execution(directory),
                                invocation_seconds=time.perf_counter() - started)
    parallel_label = f'parallel_{args.workers}'
    report = dict(execution=execution, results_equal=False, validation_passed=False)
    try:
        report['checked_trajectories_per_mode'] = compare_results(
            output / 'serial', output / parallel_label, 12, args.n_traj)
        report['results_equal'] = True
        report['elapsed_speedup'] = execution['serial']['invocation_seconds'] / execution[parallel_label]['invocation_seconds']
        report['all_paths_complete'] = all(row['incomplete_paths'] == 0 for row in execution.values())
        report['observed_parallel_execution'] = execution[parallel_label]['peak_simultaneous_cases'] > 1
        report['validation_passed'] = report['all_paths_complete'] and report['observed_parallel_execution']
    except AssertionError as error:
        report['error'] = str(error)
    (output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    if not report['validation_passed']:
        raise SystemExit('Parallel validation did not pass; inspect comparison.json.')


if __name__ == '__main__':
    main()
