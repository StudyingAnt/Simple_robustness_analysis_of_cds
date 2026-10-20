"""Run all three datasets under six enzyme profiles with one bounded array submission."""
import argparse
from copy import copy
import json
from pathlib import Path

from prepare_inputs import DATASETS, load_inputs
from prepare_profiles import digest, require
from run_experimental_apobec import run, storage_helpers


def task_cases(index, count, shards):
    require(0 < count <= shards and 0 <= index < count, 'Invalid allocation task plan')
    return [(dataset, shard) for shard in range(index, shards, count) for dataset in DATASETS]


def check_pilot(pilot_dir, inputs_dir, campaign, package_sha, seed, max_steps):
    path = Path(pilot_dir) / 'experimental_pilot_report.json'
    report = json.loads(path.read_text())
    require(report['validation_passed'] is True and report['n_trajectories'] >= 100,
            'A passing 100-path experimental pilot is required before production')
    require(report['campaign_sha256'] == digest(Path(inputs_dir) / 'inputs_manifest.json')
            and report['profiles_sha256'] == campaign['profiles_sha256']
            and report['profile_ids'] == campaign['profile_ids'], 'Pilot uses different inputs/profiles')
    for dataset in DATASETS:
        row = report['datasets'][dataset]
        require(row['results_equal'] and row['all_paths_complete'] and row['observed_parallel_execution']
                and row['audits'] == {'serial': True, 'parallel': True}, 'Experimental pilot checks did not pass')
        require(row['input_sha256'] == digest(Path(pilot_dir) / 'inputs' / f'{dataset}.fa'),
                'Selected pilot FASTA changed')
        for mode in ('serial', 'parallel'):
            configuration = json.loads((Path(pilot_dir) / dataset / mode / 'run_manifest.json').read_text())['configuration']
            require(configuration['runner_python_sha256'] == digest(Path(__file__).with_name('run_experimental_apobec.py'))
                    and configuration['profile_loader_sha256'] == digest(Path(__file__).with_name('prepare_profiles.py'))
                    and configuration['package_source_sha256'] == package_sha
                    and configuration['master_seed'] == seed and configuration['max_steps'] == max_steps,
                    'Pilot code/package/seed/cap differs from production; run a fresh pilot')


def run_task(args):
    cases = task_cases(args.task_index, args.task_count, args.shards_per_dataset)
    campaign = load_inputs(args.inputs_dir, verify_fasta=True)
    storage_helpers.load_package(args.sinabro_source)
    check_pilot(args.pilot_dir, args.inputs_dir, campaign, storage_helpers.package_digest(), args.seed, args.max_steps)
    print(f'Allocation task {args.task_index}/{args.task_count}: {cases}', flush=True)
    for dataset, shard in cases:
        case = copy(args)
        case.dataset = dataset
        case.input = Path(campaign['datasets'][dataset]['input_path'])
        case.out_dir = args.out_root / dataset
        case.profiles_dir = Path(campaign['profiles_dir'])
        case.n_shards = args.shards_per_dataset
        case.shard_index = shard
        case.expected_records = campaign['datasets'][dataset]['transcripts']
        run(case)
    print(f'Allocation task {args.task_index} complete: {len(cases)} dataset/shard cases', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-dir', type=Path, required=True)
    parser.add_argument('--pilot-dir', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--task-index', type=int, required=True)
    parser.add_argument('--task-count', type=int, default=50)
    parser.add_argument('--shards-per-dataset', type=int, default=100)
    parser.add_argument('--sinabro-source', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    run_task(parser.parse_args(argv))


if __name__ == '__main__':
    main()
