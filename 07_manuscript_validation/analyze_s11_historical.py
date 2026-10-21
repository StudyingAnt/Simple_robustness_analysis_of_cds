"""Compare recovered historical S11 estimates on the same transcript population."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import tempfile
import time

from validate_manuscript_results import digest, require, save_plot, write_csv
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
FIELDS = ['index', 'transcript_id', 'transcript_name',
          'historical_identity_rums', 'historical_blosum_rums']


def read_values(path, manifest):
    require(digest(path) == manifest['values_sha256'], 'Historical values checksum mismatch')
    raw = gzip.decompress(path.read_bytes())
    require(hashlib.sha256(raw).hexdigest() == manifest['uncompressed_csv_sha256'],
            'Historical CSV checksum mismatch')
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8'), newline=''))
    require(reader.fieldnames == FIELDS, 'Unexpected historical columns')
    rows = list(reader)
    require(len(rows) == manifest['transcripts'] and len(rows) > 0,
            'Wrong historical transcript count')
    for i, row in enumerate(rows):
        require(int(row['index']) == i, 'Historical row order changed')
        for key in FIELDS[-2:]:
            row[key] = float(row[key])
            require(np.isfinite(row[key]) and row[key] >= 0, 'Invalid historical estimate')
    require(len({r['transcript_id'] for r in rows}) == len(rows), 'Duplicate historical ID')
    require(len({r['transcript_name'] for r in rows}) == len(rows), 'Duplicate historical name')
    return rows


def check_sequences(path, rows, expected_sha):
    require(digest(path) == expected_sha, 'Human FASTA checksum mismatch')
    selected = {r['transcript_id']: r for r in rows}
    checked = {}
    with path.open() as stream:
        for record in SeqIO.parse(stream, 'fasta'):
            parts = record.id.split('|')
            if parts[0] not in selected:
                continue
            require(parts[0] not in checked, 'Duplicate selected FASTA ID')
            require(len(parts) > 4 and parts[4] == selected[parts[0]]['transcript_name'],
                    'Transcript ID/name mismatch')
            checked[parts[0]] = dict(cds_nt=len(record.seq)-2,
                sequence_sha256=hashlib.sha256(str(record.seq).encode()).hexdigest())
    require(set(checked) == set(selected), 'Missing historical transcript in human FASTA')
    return checked


def summarize(values):
    return dict(mean=float(values.mean()), median=float(np.median(values)),
                sd_across_transcripts=float(values.std(ddof=1)) if len(values) > 1 else None,
                minimum=float(values.min()), maximum=float(values.max()))


def run(args):
    started = time.perf_counter()
    manifest_sha = digest(args.manifest)
    manifest = json.loads(args.manifest.read_text())
    require(manifest['schema_version'] == 1 and manifest['historical_results_only'] is True
            and manifest['signature'] == 'SBS2', 'Unsupported historical source manifest')
    values_path = args.values or args.manifest.parent / manifest['values_filename']
    rows = read_values(values_path, manifest)
    annotation = check_sequences(args.input, rows, manifest['human_input_sha256'])
    output = args.out_dir.resolve()
    require(not output.exists(), 'Choose a new S11 output directory')
    for protected in (args.input.resolve().parent, args.manifest.resolve().parent,
                      values_path.resolve().parent):
        require(output != protected and protected not in output.parents,
                'Output cannot be nested inside an input directory')
    identity = np.array([r['historical_identity_rums'] for r in rows])
    blosum = np.array([r['historical_blosum_rums'] for r in rows])
    report = dict(analysis_passed=False, stage='s11_historical_matched_transcripts',
        historical_results_only=True, fresh_simulation=False,
        corrected_package_blosum_recalculation=False,
        n_transcripts=len(rows), input_sequences_verified=len(annotation), signature='SBS2',
        historical_script_n_sim=manifest['historical_script_n_sim'],
        historical_criterion=manifest['historical_criterion'],
        source_manifest_sha256=manifest_sha, source_values_sha256=digest(values_path),
        human_input_sha256=manifest['human_input_sha256'], script_sha256=digest(__file__),
        identity=summarize(identity), blosum=summarize(blosum),
        mean_difference_of_historical_estimates=float((blosum-identity).mean()),
        fraction_blosum_at_least_identity=float(np.mean(blosum >= identity)),
        historical_plotted_populations=manifest['historical_plotted_populations'],
        matched_comparison_populations=dict(identity=len(rows), blosum=len(rows)),
        historical_random_seeds_recovered=manifest['historical_seeds_recovered'],
        interpretation='Recovered historical simulation estimates, matched by transcript. '
            'No new trajectories, corrected-package sensitivity result, paired random paths, '
            'or measurement of functional loss is claimed.',
        runtime=dict(python=platform.python_version(), numpy=np.__version__))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.s11-historical-', dir=output.parent) as temporary:
        stage = Path(temporary)/'results'
        stage.mkdir()
        write_csv(stage/'s11_matched_transcript_values.csv', [dict(**r,
            **annotation[r['transcript_id']], historical_estimate_difference=
            r['historical_blosum_rums']-r['historical_identity_rums']) for r in rows])
        fig, ax = plt.subplots(figsize=(6, 4))
        rng = np.random.default_rng(42)  # Figure jitter only; not a simulation seed.
        for position, vector in enumerate((identity, blosum), 1):
            if len(vector) > 1 and np.ptp(vector) > 0:
                ax.violinplot([vector], positions=[position], showmedians=True)
            else:
                ax.hlines(vector[0], position-.15, position+.15)
            ax.scatter(rng.normal(position, .04, len(vector)), vector, s=2, alpha=.4)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Protein identity', 'BLOSUM80'])
        ax.set(xlabel='Non-neutral mutation criterion', ylabel='Historical mean neutral steps',
               title=f'S11: same {len(rows):,} transcripts under both criteria')
        save_plot(fig, stage, 'figureS11_historical_matched')
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(identity, blosum, s=3, alpha=.4)
        limit = max(float(identity.max()), float(blosum.max()), 1.)
        ax.plot([0, limit], [0, limit], '--', color='grey')
        ax.set(xlabel='Historical identity estimate', ylabel='Historical BLOSUM estimate',
               title='Same-transcript comparison of archived estimates')
        save_plot(fig, stage, 's11_historical_criterion_scatter')
        require(digest(args.manifest) == manifest_sha
                and digest(values_path) == report['source_values_sha256']
                and digest(args.input) == report['human_input_sha256'],
                'Historical inputs changed during analysis')
        report.update(analysis_passed=True, elapsed_seconds=time.perf_counter()-started,
            outputs=[dict(filename=p.name, sha256=digest(p)) for p in sorted(stage.iterdir())])
        (stage/'s11_historical_summary.json').write_text(
            json.dumps(report, indent=2, allow_nan=False)+'\n')
        require(not output.exists(), 'S11 output appeared during analysis')
        stage.rename(output)
    print(json.dumps({k:v for k,v in report.items() if k != 'outputs'}, indent=2), flush=True)
    print('Report:', output/'s11_historical_summary.json')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Pinned full human CDS-plus FASTA')
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path,
                        default=HERE/'source_data/s11_historical_manifest.json')
    parser.add_argument('--values', type=Path, help='Optional explicit location of the frozen gzip CSV')
    run(parser.parse_args(argv))


if __name__ == '__main__':
    main()
