"""Process the complete two-species shard plan using a fixed number of submitted jobs."""
import argparse
from copy import copy
from pathlib import Path

from bat_inputs import array_case
from extract_bat_cds import require
from run_bat_cosmic import run


def task_cases(task_index, task_count, shards_per_species):
    """Pair the species, then distribute shard pairs over fixed allocation tasks."""
    require(0 < task_count <= shards_per_species, 'Task count must be between one and shards per species')
    require(0 <= task_index < task_count, 'Task index outside the fixed submission plan')
    return [2 * shard + offset for shard in range(task_index, shards_per_species, task_count)
            for offset in (0, 1)]


def run_task(args):
    cases = task_cases(args.task_index, args.task_count, args.shards_per_species)
    plan = [dict(species=array_case(index, args.shards_per_species)[0],
                 shard=array_case(index, args.shards_per_species)[1]) for index in cases]
    print(f'Allocation task {args.task_index}/{args.task_count}: {plan}', flush=True)
    # Each case launches and joins its own worker pool before the next case begins.
    # Existing per-case provenance, seeds, shard assignment and checkpoints stay unchanged.
    for index in cases:
        case_args = copy(args)
        case_args.array_index = index
        run(case_args)
    print(f'Allocation task {args.task_index} complete: {len(cases)} species/shard cases', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--processed-root', type=Path, required=True)
    parser.add_argument('--pilot-dir', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--task-index', type=int, required=True)
    parser.add_argument('--task-count', type=int, default=50)
    parser.add_argument('--shards-per-species', type=int, default=100)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--sinabro-source', type=Path)
    run_task(parser.parse_args(argv))


if __name__ == '__main__':
    main()
