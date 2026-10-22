"""Check all six enzyme profiles on length-stratified examples from all three datasets."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_inputs import DATASETS, load_inputs
from prepare_profiles import digest, load_profiles, require

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parents[1] / '02_run_sinabro'))
sys.path.insert(0, str(HERE.parents[1] / '02_run_sinabro/author_hpc'))
from validate_parallel import summarize_execution
from validate_results import validate_run


def select_records(path, count=12):
    records = list(SeqIO.parse(path, 'fasta'))
    eligible = sorted((r for r in records if '*' not in str(r.seq[1:-1].translate())[:-1]),
                      key=lambda r: (len(r.seq) - 2, r.id))
    require(count > 0 and len(eligible) >= count, 'Insufficient eligible pilot sequences')
    return [eligible[((2 * i + 1) * len(eligible)) // (2 * count)] for i in range(count)]


def compare(left, right, n_traj, profiles, expected_records=12):
    files = {p.name for p in (left / 'transcripts').glob('*.npz')}
    require(len(files) == expected_records and files == {p.name for p in (right / 'transcripts').glob('*.npz')},
            'Pilot checkpoint sets differ')
    checked = 0
    for name in sorted(files):
        with np.load(left / 'transcripts' / name, allow_pickle=False) as a, \
             np.load(right / 'transcripts' / name, allow_pickle=False) as b:
            for key in ('accepted_steps', 'terminal_codes'):
                require(a[key].shape == b[key].shape == (n_traj, 6) and np.array_equal(a[key], b[key]),
                        f'Per-path serial/parallel mismatch: {name} {key}')
            for item in (a, b):
                require(json.loads(str(item['metadata'].item()))['signatures'] == profiles, 'Wrong enzyme profile order')
            require(json.loads(str(a['metadata'].item()))['config_sha256'] ==
                    json.loads(str(b['metadata'].item()))['config_sha256'], 'Different scientific configurations')
            checked += a['accepted_steps'].size
    return checked


def validate(args):
    require(args.workers >= 2 and args.n_traj > 0 and args.timeout > 0 and args.seed >= 0,
            'Use at least two workers, positive sizes and nonnegative seed')
    campaign = load_inputs(args.inputs_dir, verify_fasta=True)
    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'inputs').mkdir()
    report = dict(validation_passed=False, model='experimental_apobec_NNC_v1', n_trajectories=args.n_traj,
                  profile_ids=campaign['profile_ids'], profiles_sha256=campaign['profiles_sha256'],
                  campaign_sha256=digest(args.inputs_dir / 'inputs_manifest.json'), datasets={},
                  script_sha256=digest(__file__))
    try:
        for dataset in DATASETS:
            chosen = select_records(campaign['datasets'][dataset]['input_path'])
            fasta = output / 'inputs' / f'{dataset}.fa'
            SeqIO.write(chosen, fasta, 'fasta')
            target = output / dataset
            target.mkdir()
            entry = dict(input_sha256=digest(fasta), selected_cds_nt=[len(r.seq) - 2 for r in chosen],
                         selected_ids=[r.id for r in chosen], execution={}, audits={})
            report['datasets'][dataset] = entry
            for label, workers in [('serial', 1), ('parallel', args.workers)]:
                run = target / label
                command = [sys.executable, str(HERE / 'run_experimental_apobec.py'),
                    '--input', str(fasta), '--out-dir', str(run), '--profiles-dir', campaign['profiles_dir'],
                    '--dataset', dataset, '--workers', str(workers), '--n-traj', str(args.n_traj),
                    '--seed', str(args.seed), '--max-steps', str(args.max_steps), '--expected-records', '12']
                if args.sinabro_source is not None:
                    command.extend(['--sinabro-source', str(args.sinabro_source.resolve())])
                print(f'{dataset} {label}: 12 transcripts x 6 enzymes x {args.n_traj} paths', flush=True)
                started = time.perf_counter()
                subprocess.run(command, check=True, timeout=args.timeout)
                entry['execution'][label] = dict(summarize_execution(run), invocation_seconds=time.perf_counter() - started)
                audit = validate_run(run, fasta, target / f'{label}_audit', workers=1, expected_records=12)
                require(audit['validation_passed'], f'{dataset} {label} audit failed')
                entry['audits'][label] = True
            entry['checked_trajectories_per_mode'] = compare(target / 'serial', target / 'parallel', args.n_traj,
                                                            campaign['profile_ids'])
            entry['results_equal'] = True
            entry['all_paths_complete'] = all(mode['incomplete_paths'] == 0 for mode in entry['execution'].values())
            entry['observed_parallel_execution'] = entry['execution']['parallel']['peak_simultaneous_cases'] > 1
            require(entry['all_paths_complete'] and entry['observed_parallel_execution'], 'Completion/concurrency not confirmed')
        report['validation_passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        (output / 'experimental_pilot_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-dir', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--sinabro-source', type=Path)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--n-traj', type=int, default=100)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--timeout', type=int, default=3600)
    validate(parser.parse_args(argv))


if __name__ == '__main__':
    main()
