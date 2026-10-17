"""Check audited RUMS tables and draw descriptive signature distributions."""
import argparse
from contextlib import ExitStack
import csv
import hashlib
import io
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'

import numpy as np
import scipy
from scipy.stats import gaussian_kde


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verified_rows(directory, name, inventory, checked):
    """Parse exactly the CSV bytes whose hash was checked against the export."""
    require(name in inventory, f'Missing aggregation inventory entry: {name}')
    raw = (directory / name).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == inventory[name], f'Changed aggregation file: {name}')
    checked[name] = actual
    return csv.reader(io.StringIO(raw.decode('utf-8'), newline=''))


def load_matrix(directory, name, signatures, n, inventory, checked, expected_ids=None):
    reader = verified_rows(directory, name, inventory, checked)
    require(next(reader) == ['Transcript', *signatures], f'Unexpected columns: {name}')
    ids, values = [], np.empty((n, len(signatures)))
    for i, row in enumerate(reader):
        require(i < n and len(row) == len(signatures) + 1, f'Invalid row count/width: {name}')
        ids.append(row[0])
        values[i] = [float(item) for item in row[1:]]
    require(len(ids) == n and len(set(ids)) == n, f'Missing or duplicate transcript rows: {name}')
    if expected_ids is not None:
        require(ids == expected_ids, f'Transcript row correspondence mismatch: {name}')
    require(np.isfinite(values).all() and (values >= 0).all(), f'Invalid RUMS values: {name}')
    return ids, values


def signature_statistics(signatures, matrix):
    require(len(matrix) > 1, 'At least two transcripts are needed for across-transcript SD')
    averages, dispersion = matrix.mean(axis=0), matrix.std(axis=0, ddof=1)
    rows = []
    for i, name in enumerate(signatures):
        x = matrix[:, i]
        rows.append(dict(Signature=name, n_transcripts=len(x),
                         mean_rums_across_transcripts=float(averages[i]),
                         sd_rums_across_transcripts=float(dispersion[i]),
                         median_rums_across_transcripts=float(np.median(x)),
                         min_rums=float(x.min()), max_rums=float(x.max())))
    return rows


def rank_groups(rows, groups):
    lookup = {row['Signature']: row for row in rows}
    ranked = []
    for group, members in groups.items():
        for metric in ('mean', 'sd'):
            key = metric + '_rums_across_transcripts'
            ordered = sorted(members, key=lambda name: (-lookup[name][key], name))
            for rank, name in enumerate(ordered, 1):
                ranked.append(dict(group=group, ranking_metric=metric, rank=rank, **lookup[name]))
    return ranked


def density(task):
    index, values, points = task
    if values.min() == values.max():
        return index, dict(grid=np.array([values[0]]), relative_density=np.array([0.]))
    grid = np.linspace(values.min(), values.max(), points)
    estimate = gaussian_kde(values, bw_method='scott')(grid)
    return index, dict(grid=grid, relative_density=estimate / estimate.max())


def draw_panel(ax, ordered, signatures, matrix, densities, title, limit):
    for position, name in enumerate(ordered):
        index = signatures.index(name)
        x = matrix[:, index]
        curve = densities[index]
        width = .38 * curve['relative_density']
        color = '#4c78a8' if name != 'SBS2' else '#d95f02'
        if len(curve['grid']) > 1:
            ax.fill_between(curve['grid'], position-width, position+width,
                            facecolor=color, edgecolor=color, linewidth=.4, alpha=.55)
        else:
            ax.plot([x[0], x[0]], [position-.3, position+.3], color=color, linewidth=1)
        lower, median, upper = np.quantile(x, [.25, .5, .75])
        ax.plot([lower, upper], [position, position], color='#303030', linewidth=1.8)
        ax.plot(median, position, marker='o', markersize=2.4, color='white',
                markeredgecolor='#303030', markeredgewidth=.5)
    ax.set_yticks(range(len(ordered)), labels=ordered)
    ax.set_ylim(len(ordered)-.4, -.6)
    ax.set_xlim(0, limit)
    ax.set_xlabel('Transcript RUMS (mean accepted steps)', fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(axis='x', color='#dddddd', linewidth=.5)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)


def figures(stage, groups, rankings, signatures, matrix, densities, page_size, dpi):
    # Set a writable private cache before importing the noninteractive plotting backend.
    with tempfile.TemporaryDirectory(prefix='sinabro-mpl-') as cache:
        with ExitStack() as stack:
            previous = os.environ.get('MPLCONFIGDIR')
            os.environ['MPLCONFIGDIR'] = cache
            def restore_cache():
                if previous is None:
                    os.environ.pop('MPLCONFIGDIR', None)
                else:
                    os.environ['MPLCONFIGDIR'] = previous
            stack.callback(restore_cache)
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
            limit = max(float(matrix.max()) * 1.03, 0.1)
            pages = 0
            for group, members in groups.items():
                for metric in ('mean', 'sd'):
                    ordered = [r['Signature'] for r in rankings
                               if r['group'] == group and r['ranking_metric'] == metric]
                    label = 'average RUMS' if metric == 'mean' else 'across-transcript SD'
                    with PdfPages(stage / f'{group}_rank_by_{metric}.pdf') as pdf:
                        for start in range(0, len(members), page_size):
                            subset = ordered[start:start+page_size]
                            fig, ax = plt.subplots(figsize=(6, max(2.5, .30*len(subset)+1.3)))
                            draw_panel(ax, subset, signatures, matrix, densities,
                                       f'{group}: ranks {start+1}-{start+len(subset)} by {label}', limit)
                            fig.tight_layout()
                            pdf.savefig(fig)
                            fig.savefig(stage / f'{group}_{metric}_{start+1:02d}_{start+len(subset):02d}.png', dpi=dpi)
                            plt.close(fig)
                            pages += 1
            # Recreate the manuscript's top/bottom comparison from the observed 41-signature ranks.
            if 'aetiology' in groups:
                k = min(15, len(groups['aetiology']) // 2)
                if k:
                    fig, axes = plt.subplots(2, 2, figsize=(12, max(6, .6*k+2)))
                    for column, metric in enumerate(('mean', 'sd')):
                        ordered = [r['Signature'] for r in rankings
                                   if r['group'] == 'aetiology' and r['ranking_metric'] == metric]
                        for row, selected in enumerate((ordered[:k], ordered[-k:])):
                            draw_panel(axes[row, column], selected, signatures, matrix, densities,
                                       f'{"Top" if row == 0 else "Bottom"} {k} by {"mean" if metric == "mean" else "SD"}', limit)
                    fig.suptitle('Human CDS RUMS: known-aetiology signatures', fontsize=11)
                    fig.tight_layout(rect=(0, 0, 1, .97))
                    fig.savefig(stage / 'aetiology_top_bottom.pdf')
                    fig.savefig(stage / 'aetiology_top_bottom.png', dpi=dpi)
                    plt.close(fig)
            return dict(pdf_rank_pages=pages, matplotlib=matplotlib.__version__,
                        displayed_rums_range=[0, limit])


def analyze(aggregation_dir, out_dir, workers=1, page_size=15, density_points=256, dpi=200,
            expected_records=None, expected_signatures=None):
    started = time.perf_counter()
    aggregation_dir, out_dir = Path(aggregation_dir).resolve(), Path(out_dir).resolve()
    require(not out_dir.exists(), 'Output directory exists; choose a new directory')
    require(out_dir != aggregation_dir and aggregation_dir not in out_dir.parents,
            'Output must be outside the aggregation directory')
    require(workers > 0 and page_size > 0 and density_points >= 32 and dpi > 0, 'Invalid plot settings')
    manifest_path = aggregation_dir / 'aggregation_manifest.json'
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    config = manifest['configuration']
    require(manifest['aggregation_passed'] is True
            and hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest() == manifest['config_sha256'],
            'Invalid aggregation manifest/configuration')
    signatures, n = config['signatures'], manifest['transcripts']
    require(len(signatures) == len(set(signatures)) == manifest['signatures']
            and n > 1 and manifest['checkpoint_hashes_matched'] == n
            and manifest['trajectories_per_combination'] == config['n_trajectories']
            and manifest['total_trajectories'] == n * len(signatures) * config['n_trajectories'],
            'Aggregation dimensions mismatch')
    if expected_records is not None:
        require(n == expected_records, 'Unexpected transcript count')
    if expected_signatures is not None:
        require(len(signatures) == expected_signatures, 'Unexpected signature count')
    inventory = {entry['filename']: entry['sha256'] for entry in manifest['outputs']}
    require(len(inventory) == len(manifest['outputs']) and
            all(Path(name).name == name for name in inventory), 'Invalid aggregation inventory')
    checked = {}
    ids, matrix = load_matrix(aggregation_dir, 'rums_mean_by_transcript_id.csv', signatures, n, inventory, checked)
    annotation_reader = verified_rows(aggregation_dir, 'transcript_annotations.csv', inventory, checked)
    header = next(annotation_reader)
    annotations = [dict(zip(header, row)) for row in annotation_reader]
    require([r['transcript_id'] for r in annotations] == ids, 'Annotation ID correspondence mismatch')
    require(len({r['gene_id'] for r in annotations}) == manifest['unique_genes']
            and sum(r['internal_stop'] == 'True' for r in annotations)
            == manifest['internal_stop_transcripts_retained'], 'Annotation counts mismatch')
    legacy_ids = [r['transcript_name'] for r in annotations]
    groups = {'all': signatures, **manifest['signature_groups']}
    require(set(groups) <= {'all', 'noseqerr', 'aetiology'}, 'Unknown signature group')
    for group, members in groups.items():
        require(members and len(set(members)) == len(members) and set(members) <= set(signatures),
                'Unavailable or duplicate group signatures')
        name = 'all_gencode_rums_profile' + ('' if group == 'all' else '_' + group) + '.csv'
        _, legacy = load_matrix(aggregation_dir, name, members, n, inventory, checked, legacy_ids)
        require(np.array_equal(legacy, matrix[:, [signatures.index(s) for s in members]]),
                f'Canonical/legacy matrix values differ: {group}')
    del legacy
    stats = signature_statistics(signatures, matrix)
    rankings = rank_groups(stats, groups)
    source_reader = verified_rows(aggregation_dir, 'signature_summary.csv', inventory, checked)
    header = next(source_reader)
    saved_stats = {r['Signature']: r for r in (dict(zip(header, row)) for row in source_reader)}
    require(set(saved_stats) == set(signatures), 'Signature summary IDs mismatch')
    for row in stats:
        saved = saved_stats[row['Signature']]
        require(int(saved['n_transcripts']) == n, 'Summary transcript count mismatch')
        for key in row.keys() - {'Signature', 'n_transcripts'}:
            require(np.isclose(float(saved[key]), row[key], rtol=1e-10, atol=1e-12),
                    f'Signature summary mismatch: {row["Signature"]} / {key}')
    for rank in rankings:
        if rank['group'] == 'all':
            require(int(saved_stats[rank['Signature']]['rank_' + rank['ranking_metric'] + '_desc'])
                    == rank['rank'], 'Saved signature rank mismatch')
    # Hash unused uncertainty tables too: the entire published aggregation must be intact.
    for name, expected in inventory.items():
        if name not in checked:
            actual = digest_file(aggregation_dir / name)
            require(actual == expected, f'Changed aggregation file: {name}')
            checked[name] = actual
    print(f'Checked {n} transcripts x {len(signatures)} signatures; summary/rank agreement passed', flush=True)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.signature-analysis-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'analysis'
        stage.mkdir()
        with (stage / 'signature_rankings.csv').open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rankings[0]))
            writer.writeheader()
            writer.writerows(rankings)
        density_tasks = ((i, matrix[:, i], density_points) for i in range(len(signatures)))
        with ExitStack() as stack:
            pool = stack.enter_context(mp.get_context('spawn').Pool(workers)) if workers > 1 else None
            computed = pool.imap_unordered(density, density_tasks) if pool else map(density, density_tasks)
            densities = {}
            for i, curve in computed:
                densities[i] = curve
                if len(densities) % 10 == 0 or len(densities) == len(signatures):
                    print(f'Computed full-data densities {len(densities)}/{len(signatures)}', flush=True)
        plotting = figures(stage, groups, rankings, signatures, matrix, densities, page_size, dpi)
        leaders, claim_checks = {}, {}
        for group in groups:
            leaders[group] = {}
            for metric in ('mean', 'sd'):
                members = [r for r in rankings if r['group'] == group and r['ranking_metric'] == metric]
                leaders[group][metric] = dict(signature=members[0]['Signature'],
                                              value=members[0][metric + '_rums_across_transcripts'])
                sbs2 = next((r for r in members if r['Signature'] == 'SBS2'), None)
                claim_checks.setdefault(group, {})['SBS2_' + metric + '_rank'] = sbs2['rank'] if sbs2 else None
        require(manifest_path.read_bytes() == manifest_bytes, 'Aggregation manifest changed during analysis')
        result = dict(analysis_passed=True, summary_and_rank_agreement=True,
                      aggregation_directory=str(aggregation_dir),
                      aggregation_manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                      config_sha256=manifest['config_sha256'], checked_input_hashes=checked,
                      transcripts=n, unique_genes=manifest['unique_genes'], signatures=len(signatures),
                      internal_stop_transcripts_retained=manifest['internal_stop_transcripts_retained'],
                      signature_groups={g: len(m) for g, m in groups.items()}, leaders=leaders,
                      manuscript_claim_checks=claim_checks, sample_sd_ddof=1,
                      weighting='all transcript isoforms with equal transcript weight',
                      kernel='Gaussian KDE, Scott bandwidth, all transcript values; constant columns drawn as lines',
                      density_points=density_points, workers=workers, page_size=page_size, png_dpi=dpi,
                      quartiles='all observations; median point and interquartile bar',
                      subsampling=False, inferential_tests_performed=False,
                      plotting=plotting, python=platform.python_version(), numpy=np.__version__,
                      scipy=scipy.__version__, script_sha256=digest_file(__file__),
                      outputs=[dict(filename=p.name, sha256=digest_file(p), bytes=p.stat().st_size)
                               for p in sorted(stage.iterdir())],
                      elapsed_seconds=time.perf_counter()-started,
                      completed_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        (stage / 'analysis_summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        require(not out_dir.exists(), 'Output directory appeared during analysis')
        stage.rename(out_dir)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--aggregation-dir', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True, help='New directory outside the aggregation')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--page-size', type=int, default=15)
    parser.add_argument('--density-points', type=int, default=256)
    parser.add_argument('--dpi', type=int, default=200)
    parser.add_argument('--expected-records', type=int)
    parser.add_argument('--expected-signatures', type=int)
    args = parser.parse_args(argv)
    result = analyze(args.aggregation_dir, args.out_dir, args.workers, args.page_size,
                     args.density_points, args.dpi, args.expected_records, args.expected_signatures)
    print(json.dumps({key: result[key] for key in ('analysis_passed', 'summary_and_rank_agreement',
        'transcripts', 'signatures', 'leaders', 'manuscript_claim_checks', 'elapsed_seconds')}, indent=2))
    print('Report:', args.out_dir / 'analysis_summary.json')
    return 0


if __name__ == '__main__':
    mp.freeze_support()
    sys.exit(main())
