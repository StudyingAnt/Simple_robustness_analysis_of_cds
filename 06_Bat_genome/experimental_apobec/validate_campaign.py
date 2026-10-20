"""Audit every saved experimental-enzyme path in the three full datasets."""
import argparse
import json
from pathlib import Path
import sys

from prepare_inputs import DATASETS, load_inputs
from prepare_profiles import digest, require

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / '02_run_sinabro'))
from validate_results import validate_run


def validate(args):
    require(args.workers > 0 and args.n_traj > 0, 'Workers and path count must be positive')
    campaign = load_inputs(args.inputs_dir, verify_fasta=True)
    output = args.out_root.resolve()
    require(all(output != p.resolve() and p.resolve() not in output.parents
                for p in (args.inputs_dir, args.run_root, Path(campaign['profiles_dir']))),
            'Audit output must be outside simulation/input/profile directories')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(stage='audit', stage_passed=False, datasets={}, n_trajectories=args.n_traj,
                  campaign_sha256=digest(args.inputs_dir / 'inputs_manifest.json'), script_sha256=digest(__file__))
    try:
        for dataset in DATASETS:
            run = args.run_root / dataset
            config = json.loads((run / 'run_manifest.json').read_text())['configuration']
            require(config['mutation_model'] == 'experimental_apobec_NNC_v1' and config['dataset'] == dataset
                    and config['profiles_sha256'] == campaign['profiles_sha256']
                    and config['signatures'] == campaign['profile_ids']
                    and config['n_trajectories'] == args.n_traj
                    and config['criterion'] == 'protein_identity' and config['strand'] == 'both'
                    and config['strand_bias'] == .5 and config['cds_start'] == 1
                    and config['cds_end'] == 'len(sequence)-1', 'Experimental run settings mismatch')
            result = validate_run(run, campaign['datasets'][dataset]['input_path'], output / dataset,
                                  args.workers, campaign['datasets'][dataset]['transcripts'])
            require(result['input_sha256'] == campaign['datasets'][dataset]['input_sha256']
                    and result['validation_passed'], f'{dataset} full archive audit failed')
            report['datasets'][dataset] = dict(passed=True, transcripts=result['expected_transcripts'],
                paths=result['complete_paths'], report_sha256=digest(output / dataset / 'validation_summary.json'))
        report['stage_passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        (output / 'experimental_audit_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-dir', type=Path, required=True)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--workers', type=int, default=8)
    validate(parser.parse_args(argv))


if __name__ == '__main__':
    main()
