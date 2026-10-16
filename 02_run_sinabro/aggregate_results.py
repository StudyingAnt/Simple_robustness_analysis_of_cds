"""Export audited RUMS checkpoints as analysis tables without new simulations."""
import argparse
from contextlib import ExitStack
import csv
import functools
import hashlib
import io
import json
import multiprocessing as mp
from pathlib import Path
import platform
import re
import sys
import tempfile
import time

# Reuse read-only validation helpers; this also limits numerical threads to one.
from validate_results import digest_file, require
import numpy as np
import Bio
from Bio import SeqIO


METRICS = ('mean', 'sd', 'median', 'mc_se')


def read_annotations(path):
    """Keep every transcript/isoform, including records with internal stop codons."""
    rows = []
    for record in SeqIO.parse(path, 'fasta'):
        fields = record.id.split('|')
        require(len(fields) >= 6 and all(fields[i] for i in (0, 1, 4, 5)),
                'GENCODE-style FASTA identifiers with transcript/gene IDs and names are required')
        require(re.fullmatch(r'[A-Za-z0-9_.-]+', fields[0]), 'Unsafe transcript ID')
        cds = str(record.seq[1:-1])
        require(len(cds) > 0 and len(cds) % 3 == 0, 'Invalid CDS-plus length')
        rows.append(dict(transcript_id=fields[0], transcript_name=fields[4], gene_id=fields[1],
                         gene_name=fields[5], cds_nt=len(cds), fasta_id=record.id,
                         internal_stop=any(cds[i:i+3] in ('TAA', 'TAG', 'TGA')
                                           for i in range(0, len(cds)-3, 3))))
    require(rows and len({row['transcript_id'] for row in rows}) == len(rows),
            'Empty input or duplicate transcript IDs')
    require(len({row['transcript_name'] for row in rows}) == len(rows),
            'Duplicate transcript names would make the historical-format tables ambiguous')
    return sorted(rows, key=lambda row: row['transcript_id'])


def read_profile(task, configuration, fingerprint):
    path, annotation, expected_sha = task
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected_sha,
            f'Checkpoint changed since validation: {path.name}')
    with np.load(io.BytesIO(raw), allow_pickle=False) as saved:
        meta = json.loads(str(saved['metadata'].item()))
        require(meta['config_sha256'] == fingerprint and meta['fasta_id'] == annotation['fasta_id']
                and meta['transcript_id'] == annotation['transcript_id']
                and meta['signatures'] == configuration['signatures']
                and meta['n_trajectories'] == configuration['n_trajectories'],
                f'Checkpoint metadata mismatch: {path.name}')
        size, n = len(configuration['signatures']), configuration['n_trajectories']
        require(np.array_equal(saved['n_complete'], np.full(size, n))
                and np.array_equal(saved['n_incomplete'], np.zeros(size)),
                f'Incomplete checkpoint: {path.name}')
        result = {name: saved[name] for name in METRICS}
        for name, values in result.items():
            require(values.shape == (size,), f'Invalid {name} shape: {path.name}')
            if name in ('sd', 'mc_se') and n == 1:
                require(np.isnan(values).all(), f'Expected undefined {name} for one path')
            else:
                require(np.isfinite(values).all() and (values >= 0).all(),
                        f'Invalid {name} values: {path.name}')
    return result


def summarize_signatures(signatures, means, mc_se):
    """Transcript-level summaries; path dispersion and biological replication differ."""
    n = len(means)
    averages = means.mean(axis=0)
    dispersion = means.std(axis=0, ddof=1) if n > 1 else np.full(len(signatures), np.nan)
    rank_mean = {i: rank for rank, i in enumerate(
        sorted(range(len(signatures)), key=lambda i: (-averages[i], signatures[i])), 1)}
    rank_sd = {i: rank for rank, i in enumerate(
        sorted(range(len(signatures)), key=lambda i: (-dispersion[i] if n > 1 else 0, signatures[i])), 1)}
    return [dict(Signature=name, n_transcripts=n, mean_rums_across_transcripts=averages[i],
                 sd_rums_across_transcripts=dispersion[i],
                 median_rums_across_transcripts=np.median(means[:, i]),
                 min_rums=means[:, i].min(), max_rums=means[:, i].max(),
                 median_mc_se_of_transcript_mean=np.median(mc_se[:, i]),
                 rank_mean_desc=rank_mean[i], rank_sd_desc=rank_sd[i] if n > 1 else '')
            for i, name in enumerate(signatures)]


def aggregate(run_dir, validation_dir, input_path, out_dir, workers=1,
              signature_groups=None, expected_records=None):
    started = time.perf_counter()
    run_dir, validation_dir, out_dir = map(lambda p: Path(p).resolve(),
                                          (run_dir, validation_dir, out_dir))
    require(workers > 0, 'Worker count must be positive')
    require(not out_dir.exists(), 'Output directory already exists; choose a new directory')
    require(all(out_dir != p and p not in out_dir.parents for p in (run_dir, validation_dir)),
            'Output directory must be outside the simulation and validation directories')
    manifest_bytes = (run_dir / 'run_manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    config, fingerprint = manifest['configuration'], manifest['config_sha256']
    require(hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest() == fingerprint,
            'Invalid simulation configuration fingerprint')
    report_bytes = (validation_dir / 'validation_summary.json').read_bytes()
    report = json.loads(report_bytes)
    require(report['validation_passed'] is True and report['all_paths_complete'] is True
            and report['issues'] == [], 'A passing, all-complete full audit is required')
    require(report['config_sha256'] == fingerprint and report['configuration'] == config
            and report['manifest_sha256'] == hashlib.sha256(manifest_bytes).hexdigest(),
            'Validation report belongs to a different simulation manifest')
    input_path = Path(input_path or report['input_path']).resolve()
    require(digest_file(input_path) == config['input_sha256'] == report['input_sha256'],
            'Input FASTA hash mismatch')
    annotations = read_annotations(input_path)
    ids = {row['transcript_id'] for row in annotations}
    n, signatures = config['n_trajectories'], config['signatures']
    require(config['runner'] == 'rums-v3' and n > 0 and len(signatures) > 0
            and len(set(signatures)) == len(signatures), 'Unsupported simulation configuration')
    if expected_records is not None:
        require(len(annotations) == expected_records, 'Unexpected transcript count')
    require(report['expected_transcripts'] == report['structurally_valid_files']
            == report['observed_transcript_files'] == len(annotations)
            and report['completed_shards'] == report['expected_shards'] == config['shards'],
            'Audit dimensions differ from the simulation')
    require(report['expected_paths'] == report['checked_paths'] == report['complete_paths']
            == len(annotations) * n * len(signatures), 'Audit path count mismatch')
    require(all(report[key] == 0 for key in ('missing_transcripts', 'extra_transcripts',
                'invalid_files', 'no_available_event_paths', 'max_steps_paths')),
            'Audit contains missing, invalid or incomplete results')
    require({p.name for p in (run_dir / 'transcripts').glob('*.npz')}
            == {identifier + '.npz' for identifier in ids}, 'Checkpoint file set changed since audit')
    inventory_path = validation_dir / 'transcript_audit.csv'
    inventory_sha = digest_file(inventory_path)
    hashes = {}
    with inventory_path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            identifier = row['transcript_id']
            require(identifier not in hashes and row['structurally_valid'] == 'True'
                    and row['error'] == '' and int(row['complete_paths']) == n * len(signatures)
                    and int(row['no_available_event_paths']) == int(row['max_steps_paths']) == 0
                    and re.fullmatch(r'[0-9a-f]{64}', row['sha256']), 'Invalid audit inventory row')
            hashes[identifier] = row['sha256']
    require(set(hashes) == ids, 'Audit inventory does not cover the input exactly')
    groups, groups_sha = {}, None
    if signature_groups is not None:
        signature_groups = Path(signature_groups).resolve()
        groups_sha = digest_file(signature_groups)
        groups = json.loads(signature_groups.read_text())['groups']
        require(set(groups) <= {'noseqerr', 'aetiology'}, 'Unknown manuscript signature group')
        for name, members in groups.items():
            require(members and len(set(members)) == len(members) and set(members) <= set(signatures),
                    f'Unavailable or duplicate signatures in group: {name}')
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.rums-aggregation-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'tables'
        stage.mkdir()
        means = np.empty((len(annotations), len(signatures)))
        errors = np.empty_like(means)
        with ExitStack() as stack:
            def writer(filename, header):
                stream = stack.enter_context((stage / filename).open('x', newline=''))
                result = csv.writer(stream)
                result.writerow(header)
                return result
            canonical = {metric: writer(f'rums_{metric}_by_transcript_id.csv', ['Transcript', *signatures])
                         for metric in METRICS}
            legacy = writer('all_gencode_rums_profile.csv', ['Transcript', *signatures])
            subsets = {name: writer(f'all_gencode_rums_profile_{name}.csv', ['Transcript', *members])
                       for name, members in groups.items()}
            indices = {name: [signatures.index(member) for member in members]
                       for name, members in groups.items()}
            annotation_writer = writer('transcript_annotations.csv', list(annotations[0]))
            pool = stack.enter_context(mp.get_context('spawn').Pool(workers)) if workers > 1 else None
            tasks = ((run_dir / 'transcripts' / (a['transcript_id'] + '.npz'), a, hashes[a['transcript_id']])
                     for a in annotations)
            load = functools.partial(read_profile, configuration=config, fingerprint=fingerprint)
            results = pool.imap(load, tasks, chunksize=8) if pool else map(load, tasks)
            for index, (annotation, values) in enumerate(zip(annotations, results)):
                for metric, output in canonical.items():
                    output.writerow([annotation['transcript_id'], *values[metric]])
                legacy.writerow([annotation['transcript_name'], *values['mean']])
                for name, output in subsets.items():
                    output.writerow([annotation['transcript_name'], *values['mean'][indices[name]]])
                annotation_writer.writerow(list(annotation.values()))
                means[index], errors[index] = values['mean'], values['mc_se']
                if (index + 1) % 1000 == 0 or index + 1 == len(annotations):
                    print(f'Exported {index+1}/{len(annotations)} audited transcripts', flush=True)
            summary_rows = summarize_signatures(signatures, means, errors)
            output = writer('signature_summary.csv', list(summary_rows[0]))
            for row in summary_rows:
                output.writerow(list(row.values()))
        require((run_dir / 'run_manifest.json').read_bytes() == manifest_bytes
                and (validation_dir / 'validation_summary.json').read_bytes() == report_bytes
                and digest_file(inventory_path) == inventory_sha, 'Inputs changed during export')
        outputs = [dict(filename=p.name, bytes=p.stat().st_size, sha256=digest_file(p))
                   for p in sorted(stage.iterdir())]
        result = dict(aggregation_passed=True, config_sha256=fingerprint, configuration=config,
                      run_directory=str(run_dir), validation_directory=str(validation_dir),
                      input_path=str(input_path), input_sha256=config['input_sha256'],
                      validation_report_sha256=hashlib.sha256(report_bytes).hexdigest(),
                      audit_inventory_sha256=inventory_sha, checkpoint_hashes_matched=len(annotations),
                      script_sha256=digest_file(__file__), workers=workers,
                      python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__,
                      transcripts=len(annotations), unique_genes=len({r['gene_id'] for r in annotations}),
                      internal_stop_transcripts_retained=sum(r['internal_stop'] for r in annotations),
                      signatures=len(signatures), trajectories_per_combination=n,
                      total_trajectories=len(annotations) * n * len(signatures),
                      signature_groups=groups, signature_groups_sha256=groups_sha,
                      row_order='transcript_id ascending',
                      legacy_Transcript_column='GENCODE transcript_name (FASTA field 5)',
                      canonical_Transcript_column='versioned GENCODE transcript_id (FASTA field 1)',
                      weighting='all transcript isoforms, equal transcript weight; no gene-level collapse',
                      path_sd_ddof=1, across_transcript_sd_ddof=1, outputs=outputs,
                      elapsed_seconds=time.perf_counter()-started,
                      completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        (stage / 'aggregation_manifest.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        require(not out_dir.exists(), 'Output directory appeared during export')
        stage.rename(out_dir)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--validation-dir', type=Path, required=True, help='Passing full audit directory')
    parser.add_argument('--input', type=Path, help='Default: FASTA path in the audit report')
    parser.add_argument('--out-dir', type=Path, required=True, help='New directory for analysis tables')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--signature-groups', type=Path, help='Optional manuscript subset JSON')
    parser.add_argument('--expected-records', type=int)
    args = parser.parse_args(argv)
    if args.workers < 1 or (args.expected_records is not None and args.expected_records < 1):
        parser.error('Worker and expected-record counts must be positive')
    result = aggregate(args.run_dir, args.validation_dir, args.input, args.out_dir,
                       args.workers, args.signature_groups, args.expected_records)
    print(json.dumps({key: result[key] for key in ('aggregation_passed', 'transcripts', 'unique_genes',
        'internal_stop_transcripts_retained', 'signatures', 'trajectories_per_combination',
        'total_trajectories', 'checkpoint_hashes_matched', 'elapsed_seconds')}, indent=2))
    print('Manifest:', args.out_dir / 'aggregation_manifest.json')
    return 0


if __name__ == '__main__':
    mp.freeze_support()
    sys.exit(main())
