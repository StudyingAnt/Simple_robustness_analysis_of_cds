"""Export audited bat checkpoints without inventing gene identifiers."""
import argparse
from contextlib import ExitStack
import csv
import functools
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import re
import sys
import tempfile
import time

from bat_inputs import HERE, SPECIES, check_config, read_inputs
from extract_bat_cds import digest, require

sys.path.insert(0, str(HERE.parent / '02_run_sinabro'))
from aggregate_results import METRICS, read_profile, summarize_signatures
import numpy as np


def aggregate_bat(run_dir, validation_dir, processed_dir, out_dir, species,
                  workers=1, n_traj=1000, input_path=None, allow_subset=False,
                  cohort_manifest=HERE / 'cohorts/manifest.json', configuration_checker=check_config):
    started = time.perf_counter()
    run_dir, validation_dir, processed_dir, out_dir = [Path(p).resolve() for p in
                                                    (run_dir, validation_dir, processed_dir, out_dir)]
    require(workers > 0, 'Workers must be positive')
    require(not out_dir.exists(), 'Output directory already exists')
    require(all(out_dir != p and p not in out_dir.parents for p in (run_dir, validation_dir, processed_dir)),
            'Output must be outside run, audit and processed input directories')
    inputs = read_inputs(processed_dir, species, input_path, allow_subset, cohort_manifest)
    extraction_sha = digest(inputs['manifest_path'])
    annotations = inputs['annotations']
    ids = {row['transcript_id'] for row in annotations}
    manifest_path = run_dir / 'run_manifest.json'
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    config, fingerprint = manifest['configuration'], manifest['config_sha256']
    require(hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest() == fingerprint,
            'Invalid run configuration fingerprint')
    configuration_checker(config, inputs, n_traj)
    report_path = validation_dir / 'validation_summary.json'
    report_bytes = report_path.read_bytes()
    report = json.loads(report_bytes)
    require(report['validation_passed'] is True and report['all_paths_complete'] is True and report['issues'] == [],
            'A passing complete archive audit is required')
    require(report['configuration'] == config and report['config_sha256'] == fingerprint
            and report['manifest_sha256'] == hashlib.sha256(manifest_bytes).hexdigest()
            and report['input_sha256'] == inputs['input_sha256'], 'Audit belongs to another run/input')
    require(report['expected_transcripts'] == report['observed_transcript_files']
            == report['structurally_valid_files'] == len(ids)
            and report['completed_shards'] == report['expected_shards'] == config['shards'],
            'Audit dimensions differ from the input/run')
    total_paths = len(ids) * n_traj * len(config['signatures'])
    require(report['expected_paths'] == report['checked_paths'] == report['complete_paths'] == total_paths,
            'Audit trajectory counts differ')
    require(all(report[key] == 0 for key in ('missing_transcripts', 'extra_transcripts', 'invalid_files',
                                           'no_available_event_paths', 'max_steps_paths')),
            'Audit contains invalid/missing/incomplete results')
    require({p.name for p in (run_dir / 'transcripts').glob('*.npz')} == {name + '.npz' for name in ids},
            'Checkpoint set changed since audit')
    inventory = validation_dir / 'transcript_audit.csv'
    inventory_sha = digest(inventory)
    hashes = {}
    with inventory.open(newline='') as stream:
        for row in csv.DictReader(stream):
            name = row['transcript_id']
            require(name not in hashes and row['structurally_valid'] == 'True' and row['error'] == ''
                    and int(row['complete_paths']) == n_traj * len(config['signatures'])
                    and int(row['no_available_event_paths']) == int(row['max_steps_paths']) == 0
                    and re.fullmatch('[0-9a-f]{64}', row['sha256']), 'Invalid audit inventory entry')
            hashes[name] = row['sha256']
    require(set(hashes) == ids, 'Audit inventory does not cover the input exactly')
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.bat-export-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'tables'
        stage.mkdir()
        means = np.empty((len(annotations), len(config['signatures'])))
        mc_se = np.empty_like(means)
        with ExitStack() as stack:
            def writer(name, header):
                stream = stack.enter_context((stage / name).open('w', newline=''))
                result = csv.writer(stream)
                result.writerow(header)
                return result
            outputs = {metric: writer(f'rums_{metric}_by_transcript_id.csv', ['Transcript', *config['signatures']])
                       for metric in METRICS}
            legacy = writer(f'{species}_rums_profile.csv', ['Transcript', *config['signatures']])
            annotation_writer = writer('transcript_annotations.csv', list(annotations[0]))
            pool = stack.enter_context(mp.get_context('spawn').Pool(workers)) if workers > 1 else None
            # The generic read-only checkpoint loader needs only IDs and the bound FASTA ID.
            tasks = ((run_dir / 'transcripts' / (row['transcript_id'] + '.npz'),
                      dict(transcript_id=row['transcript_id'], fasta_id=row['transcript_id']), hashes[row['transcript_id']])
                     for row in annotations)
            reader = functools.partial(read_profile, configuration=config, fingerprint=fingerprint)
            results = pool.imap(reader, tasks, chunksize=8) if pool else map(reader, tasks)
            for index, (annotation, values) in enumerate(zip(annotations, results)):
                for metric, output in outputs.items():
                    output.writerow([annotation['transcript_id'], *values[metric]])
                legacy.writerow([annotation['transcript_id'], *values['mean']])
                annotation_writer.writerow(list(annotation.values()))
                means[index], mc_se[index] = values['mean'], values['mc_se']
                if (index + 1) % 1000 == 0 or index + 1 == len(annotations):
                    print(f'{species}: exported {index+1}/{len(annotations)} audited transcripts', flush=True)
            summary = summarize_signatures(config['signatures'], means, mc_se)
            output = writer('signature_summary.csv', list(summary[0]))
            output.writerows(list(row.values()) for row in summary)
        require(manifest_path.read_bytes() == manifest_bytes and report_path.read_bytes() == report_bytes
                and digest(inventory) == inventory_sha and digest(inputs['manifest_path']) == extraction_sha
                and digest(inputs['input_path']) == inputs['input_sha256'], 'Inputs changed during export')
        files = [dict(filename=p.name, bytes=p.stat().st_size, sha256=digest(p)) for p in sorted(stage.iterdir())]
        result = dict(aggregation_passed=True, species=species, population=inputs['scope'],
                      full_historical_transcripts=inputs['full_transcripts'], transcripts=len(annotations),
                      gene_mapping_available=False, weighting='all retained transcript isoforms with equal transcript weight',
                      internal_stop_transcripts_retained=sum(int(r['internal_stop_count']) > 0 for r in annotations),
                      signatures=len(config['signatures']), trajectories_per_combination=n_traj,
                      total_trajectories=total_paths, checkpoint_hashes_matched=len(annotations),
                      config_sha256=fingerprint, configuration=config,
                      input_sha256=inputs['input_sha256'], extraction_manifest_sha256=extraction_sha,
                      audit_report_sha256=hashlib.sha256(report_bytes).hexdigest(), audit_inventory_sha256=inventory_sha,
                      run_directory=str(run_dir), validation_directory=str(validation_dir),
                      input_path=str(inputs['input_path']), script_sha256=digest(__file__),
                      input_checker_sha256=digest(HERE / 'bat_inputs.py'),
                      checkpoint_reader_sha256=digest(HERE.parent / '02_run_sinabro/aggregate_results.py'),
                      path_sd_ddof=1, across_transcript_sd_ddof=1, outputs=files,
                      elapsed_seconds=time.perf_counter() - started)
        (stage / 'aggregation_manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        require(not out_dir.exists(), 'Output appeared during export')
        stage.rename(out_dir)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--validation-dir', type=Path, required=True)
    parser.add_argument('--processed-dir', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--species', choices=SPECIES, required=True)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--input', type=Path, help='Only with --allow-subset; default full validated cohort')
    parser.add_argument('--allow-subset', action='store_true', help='Explicitly label/export a bounded input subset')
    args = parser.parse_args(argv)
    result = aggregate_bat(args.run_dir, args.validation_dir, args.processed_dir, args.out_dir, args.species,
                           args.workers, args.n_traj, args.input, args.allow_subset)
    print(json.dumps({key: result[key] for key in ('aggregation_passed', 'species', 'population', 'transcripts',
        'internal_stop_transcripts_retained', 'signatures', 'total_trajectories')}, indent=2))


if __name__ == '__main__':
    mp.freeze_support()
    main()
