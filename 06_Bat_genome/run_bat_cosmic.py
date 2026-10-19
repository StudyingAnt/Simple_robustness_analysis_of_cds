"""Run one species/shard of the frozen-cohort, eight-signature bat COSMIC analysis."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from bat_inputs import HERE, SIGNATURES, array_case, check_pilot, read_inputs
from extract_bat_cds import digest, require

RUNNER_DIR = HERE.parent / '02_run_sinabro'
sys.path.insert(0, str(RUNNER_DIR))
from run_rums import atomic_json, exclusive_lock


def run(args):
    require(args.workers > 0 and args.n_traj > 0 and args.seed >= 0 and args.max_steps > 0,
            'Invalid workers/paths/seed/cap')
    species, shard = array_case(args.array_index, args.shards_per_species)
    processed = args.processed_root.resolve()
    pilot_path = check_pilot(args.pilot_dir, processed)
    inputs = read_inputs(processed / species, species)
    output = args.out_root.resolve() / species
    require(output != processed and processed not in output.parents, 'Run output must be outside processed inputs')
    require(output != args.pilot_dir.resolve() and args.pilot_dir.resolve() not in output.parents,
            'Run output must be outside pilot outputs')
    output.mkdir(parents=True, exist_ok=True)
    provenance = dict(species=species, population=inputs['scope'],
                      expected_transcripts=len(inputs['annotations']), input_sha256=inputs['input_sha256'],
                      extraction_manifest_sha256=digest(inputs['manifest_path']),
                      cohort_sha256=inputs['manifest']['cohort_sha256'],
                      pilot_report_sha256=digest(pilot_path),
                      driver_sha256=digest(__file__), input_checker_sha256=digest(HERE / 'bat_inputs.py'),
                      signatures=SIGNATURES, n_trajectories=args.n_traj, seed=args.seed,
                      max_steps=args.max_steps, shards=args.shards_per_species)
    with exclusive_lock(output / '.bat-input.lock'):
        path = output / 'bat_input_provenance.json'
        if path.exists():
            require(json.loads(path.read_text()) == provenance,
                    'Bat output provenance/configuration changed; choose a new run root')
        else:
            require(not (output / 'run_manifest.json').exists(), 'Existing run has no bat production provenance')
            atomic_json(path, provenance)
    command = [sys.executable, str(RUNNER_DIR / 'run_rums.py'),
               '--input', str(inputs['input_path']), '--out-dir', str(output),
               '--n-traj', str(args.n_traj), '--seed', str(args.seed), '--max-steps', str(args.max_steps),
               '--workers', str(args.workers), '--n-shards', str(args.shards_per_species),
               '--shard-index', str(shard), '--expected-records', str(len(inputs['annotations'])),
               '--start-method', 'spawn', '--signatures', *SIGNATURES]
    if args.sinabro_source is not None:
        command.extend(['--sinabro-source', str(args.sinabro_source.resolve())])
    print(f'Array {args.array_index}: {species}, shard {shard}/{args.shards_per_species}, '
          f'{len(inputs["annotations"])} transcripts, {args.n_traj} paths per combination', flush=True)
    subprocess.run(command, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--processed-root', type=Path, required=True)
    parser.add_argument('--pilot-dir', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--array-index', type=int, required=True)
    parser.add_argument('--shards-per-species', type=int, default=100)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--sinabro-source', type=Path)
    run(parser.parse_args(argv))


if __name__ == '__main__':
    main()
