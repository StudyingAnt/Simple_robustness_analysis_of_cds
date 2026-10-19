"""Audit or export both species of a completed full-cohort bat COSMIC run."""
import argparse
import json
from pathlib import Path
import sys

from bat_inputs import HERE, SPECIES, check_config, read_inputs
from extract_bat_cds import digest, require
from aggregate_bat_results import aggregate_bat

sys.path.insert(0, str(HERE.parent / '02_run_sinabro'))
from validate_results import validate_run


def process(args):
    require(args.workers > 0 and args.n_traj > 0, 'Workers and paths must be positive')
    processed, runs, output = [p.resolve() for p in (args.processed_root, args.run_root, args.out_root)]
    protected = [processed, runs]
    if args.audit_root is not None:
        protected.append(args.audit_root.resolve())
    require(all(output != p and p not in output.parents for p in protected), 'Output must be outside inputs')
    require(args.stage != 'export' or args.audit_root is not None, 'Export requires --audit-root')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(stage=args.stage, stage_passed=False, n_trajectories=args.n_traj, species={},
                  run_root=str(runs), processed_root=str(processed), script_sha256=digest(__file__))
    try:
        for species in SPECIES:
            inputs = read_inputs(processed / species, species)
            config = json.loads((runs / species / 'run_manifest.json').read_text())['configuration']
            check_config(config, inputs, args.n_traj)
            if args.stage == 'audit':
                result = validate_run(runs / species, inputs['input_path'], output / species,
                                      args.workers, len(inputs['annotations']))
                require(result['validation_passed'], f'{species} archive audit failed; inspect {output / species}')
                report['species'][species] = dict(passed=True, paths=result['complete_paths'],
                    transcripts=result['expected_transcripts'], report_sha256=digest(output / species / 'validation_summary.json'))
            else:
                result = aggregate_bat(runs / species, args.audit_root / species, processed / species,
                                       output / species, species, args.workers, args.n_traj)
                report['species'][species] = dict(passed=True, paths=result['total_trajectories'],
                    transcripts=result['transcripts'], report_sha256=digest(output / species / 'aggregation_manifest.json'))
        report['stage_passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        (output / f'{args.stage}_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('audit', 'export'), required=True)
    parser.add_argument('--processed-root', type=Path, required=True)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--audit-root', type=Path)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--n-traj', type=int, default=1000)
    process(parser.parse_args(argv))


if __name__ == '__main__':
    main()
