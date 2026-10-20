"""Export the three audited experimental APOBEC datasets without resimulation."""
import argparse
import csv
import functools
import json
from pathlib import Path
import sys
import tempfile
import time

from Bio import SeqIO

from prepare_inputs import DATASETS, load_inputs
from prepare_profiles import digest, load_profiles, require

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1] / '02_run_sinabro'))
from aggregate_results import aggregate
from aggregate_bat_results import aggregate_bat
from bat_inputs import read_inputs


def check_configuration(config, inputs, n_traj, campaign, dataset):
    require(config['runner'] == 'rums-v3' and config['mutation_model'] == 'experimental_apobec_NNC_v1'
            and config['dataset'] == dataset and config['input_sha256'] == inputs['input_sha256'],
            'Experimental run belongs to a different model/dataset/input')
    require(config['signatures'] == campaign['profile_ids']
            and config['profiles_sha256'] == campaign['profiles_sha256']
            and config['profile_preparation_sha256'] == campaign['profile_manifest_sha256'],
            'Experimental profile identity/order/hash mismatch')
    require(config['n_trajectories'] == n_traj and n_traj > 0
            and config['criterion'] == 'protein_identity' and config['strand'] == 'both'
            and config['strand_bias'] == .5 and config['cds_start'] == 1
            and config['cds_end'] == 'len(sequence)-1' and config['translation_table'] == 1
            and config['internal_stops_retained_as_preprocessed'] is True
            and config['cosmic_version'] is None and config['genome_ref'] is None,
            'Unexpected experimental scientific configuration')


def export_input(campaign, dataset, manifest, allow_subset):
    full = Path(campaign['datasets'][dataset]['input_path'])
    actual = Path(manifest['input_path']) if allow_subset else full
    if dataset != 'human':
        processed = Path(campaign['datasets'][dataset]['extraction_manifest_path']).parent
        return read_inputs(processed, dataset, input_path=actual, allow_subset=allow_subset)
    scope = 'full_historical_cohort'
    if actual.resolve() != full.resolve():
        require(allow_subset, 'A bounded export requires explicit --allow-subset')
        records = list(SeqIO.parse(actual, 'fasta'))
        reference = {r.id: r for r in SeqIO.parse(full, 'fasta')}
        require(records and len(records) == len({r.id for r in records})
                and all(r.id in reference and r.description == reference[r.id].description
                        and r.seq == reference[r.id].seq for r in records),
                'Human subset is not identical to the frozen input sequences')
        scope = 'explicit_subset_of_historical_cohort'
    return dict(input_path=actual, input_sha256=digest(actual), scope=scope,
                transcripts=sum(1 for _ in SeqIO.parse(actual, 'fasta')))


def export(args):
    require(args.workers > 0 and args.n_traj > 0, 'Workers and trajectory count must be positive')
    campaign = load_inputs(args.inputs_dir, verify_fasta=True)
    profiles = load_profiles(campaign['profiles_dir'])
    profile_by_id = {row['id']: row for row in profiles['profiles']}
    runs, audits, output = [p.resolve() for p in (args.run_root, args.audit_root, args.out_root)]
    protected = [runs, audits, args.inputs_dir.resolve(), Path(campaign['profiles_dir']).resolve()]
    protected.extend(Path(row['input_path']).resolve().parent for row in campaign['datasets'].values())
    require(not output.exists(), 'Output directory already exists; choose a new directory')
    require(all(output != p and p not in output.parents for p in protected), 'Output must be outside source directories')
    campaign_path = args.inputs_dir / 'inputs_manifest.json'
    campaign_sha = digest(campaign_path)
    audit_path = audits / 'experimental_audit_manifest.json'
    audit_bytes = audit_path.read_bytes()
    audit = json.loads(audit_bytes)
    require(audit['stage'] == 'audit' and audit['stage_passed'] is True
            and audit['n_trajectories'] == args.n_traj and set(audit['datasets']) == set(DATASETS)
            and audit['campaign_sha256'] == campaign_sha, 'A passing audit for this campaign/path count is required')
    started = time.perf_counter()
    result = dict(stage='export', stage_passed=False, n_trajectories=args.n_traj, datasets={},
                  campaign_sha256=campaign_sha, audit_manifest_sha256=digest(audit_path),
                  profiles_sha256=campaign['profiles_sha256'], script_sha256=digest(__file__),
                  run_root=str(runs), audit_root=str(audits), bounded_subset_requested=args.allow_subset)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.experimental-apobec-export-', dir=output.parent) as temporary:
        stage = Path(temporary) / 'tables'
        stage.mkdir()
        summary_rows = []
        common_settings = None
        for dataset in DATASETS:
            row = audit['datasets'][dataset]
            require(row['passed'] is True and row['report_sha256'] == digest(audits / dataset / 'validation_summary.json'),
                    f'{dataset} audit summary hash/status mismatch')
            manifest = json.loads((runs / dataset / 'run_manifest.json').read_text())
            config = manifest['configuration']
            inputs = export_input(campaign, dataset, manifest, args.allow_subset)
            check_configuration(config, inputs, args.n_traj, campaign, dataset)
            settings = {key: config[key] for key in ('master_seed', 'max_steps', 'runner_python_sha256',
                'package_source_sha256', 'profile_loader_sha256', 'python', 'numpy', 'biopython', 'sinabro')}
            require(common_settings is None or settings == common_settings,
                    'Datasets use different simulation source/package/environment/seed/cap')
            common_settings = settings
            count = inputs['transcripts'] if dataset == 'human' else len(inputs['annotations'])
            require(row['transcripts'] == count and row['paths'] == count * args.n_traj * 6,
                    f'{dataset} campaign audit dimensions mismatch')
            if not args.allow_subset:
                require(count == campaign['datasets'][dataset]['transcripts'], 'Incomplete production population')
            print(f'{dataset}: exporting {count} audited transcripts under all six enzyme profiles', flush=True)
            if dataset == 'human':
                exported = aggregate(runs / dataset, audits / dataset, inputs['input_path'], stage / dataset,
                                     workers=args.workers, expected_records=count)
            else:
                processed = Path(campaign['datasets'][dataset]['extraction_manifest_path']).parent
                checker = functools.partial(check_configuration, campaign=campaign, dataset=dataset)
                exported = aggregate_bat(runs / dataset, audits / dataset, processed, stage / dataset, dataset,
                    workers=args.workers, n_traj=args.n_traj, input_path=inputs['input_path'],
                    allow_subset=args.allow_subset, configuration_checker=checker)
            result['datasets'][dataset] = dict(passed=True, transcripts=count, paths=exported['total_trajectories'],
                population=inputs['scope'], internal_stop_transcripts_retained=exported['internal_stop_transcripts_retained'],
                report_sha256=digest(stage / dataset / 'aggregation_manifest.json'))
            with (stage / dataset / 'signature_summary.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            require([row['Signature'] for row in rows] == campaign['profile_ids'], 'Export changed enzyme order')
            for row in rows:
                profile = profile_by_id[row['Signature']]
                relation = 'human_control' if dataset == 'human' else (
                    'same_species' if dataset == profile['enzyme_species'] else 'other_species')
                summary_rows.append(dict(dataset=dataset, **row, enzyme_species=profile['enzyme_species'],
                                         sample_code=profile['sample_code'], relation=relation))
        require(digest(campaign_path) == campaign_sha and audit_path.read_bytes() == audit_bytes
                and all(digest(audits / dataset / 'validation_summary.json') == audit['datasets'][dataset]['report_sha256']
                        for dataset in DATASETS), 'Campaign/audit changed during export')
        load_profiles(campaign['profiles_dir'])
        with (stage / 'profile_comparison_summary.csv').open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
            writer.writeheader()
            writer.writerows(summary_rows)
        result.update(stage_passed=True, total_transcripts=sum(row['transcripts'] for row in result['datasets'].values()),
                      total_trajectories=sum(row['paths'] for row in result['datasets'].values()),
                      summary_rows=len(summary_rows), elapsed_seconds=time.perf_counter() - started,
                      weighting='all retained transcript isoforms, equal transcript weight; no gene/orthologue collapse',
                      outputs=[dict(filename=str(p.relative_to(stage)), bytes=p.stat().st_size, sha256=digest(p))
                               for p in sorted(stage.rglob('*')) if p.is_file()])
        (stage / 'export_manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        require(not output.exists(), 'Output appeared during export')
        stage.rename(output)
    print(json.dumps({key: value for key, value in result.items() if key != 'outputs'}, indent=2), flush=True)
    print('Manifest:', output / 'export_manifest.json', flush=True)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-dir', type=Path, required=True)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--audit-root', type=Path, required=True)
    parser.add_argument('--out-root', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--allow-subset', action='store_true', help='Explicit bounded demonstration; never the default full export')
    export(parser.parse_args(argv))


if __name__ == '__main__':
    main()
