"""Recreate original-paper downstream analyses from audited exported tables."""
import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import tempfile
import time

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy import stats
import scipy
from Bio import SeqIO
if 'MPLCONFIGDIR' not in os.environ:
    _plot_cache = tempfile.TemporaryDirectory(prefix='sinabro-plot-cache-')
    os.environ['MPLCONFIGDIR'] = _plot_cache.name
    os.environ.setdefault('XDG_CACHE_HOME', _plot_cache.name)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SPECIES = ('human', 'molMol2', 'myoMyo6')
COUNTS = (73214, 31185, 23213)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_csv(path, rows):
    rows = list(rows)
    require(rows, f'No rows for {path.name}')
    with Path(path).open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_export(folder):
    folder = Path(folder).resolve()
    manifest = json.loads((folder / 'aggregation_manifest.json').read_text())
    require(manifest['aggregation_passed'] is True, f'Unsuccessful export: {folder}')
    fingerprint = hashlib.sha256(json.dumps(manifest['configuration'], sort_keys=True).encode()).hexdigest()
    require(manifest['config_sha256'] == fingerprint, 'Export configuration fingerprint mismatch')
    for entry in manifest['outputs']:
        path = folder / entry['filename']
        require(path.resolve().parent == folder and digest(path) == entry['sha256'], f'Changed export: {path}')
    signatures = manifest['configuration']['signatures']
    with (folder / 'rums_mean_by_transcript_id.csv').open(newline='') as stream:
        reader = csv.reader(stream)
        require(next(reader) == ['Transcript', *signatures], 'Mean table header/order mismatch')
        rows = list(reader)
    ids = [row[0] for row in rows]
    matrix = np.array([[float(x) for x in row[1:]] for row in rows])
    require(len(ids) == len(set(ids)) == manifest['transcripts']
            and matrix.shape == (len(ids), len(signatures)) and np.isfinite(matrix).all()
            and (matrix >= 0).all(), 'Invalid transcript matrix')
    require(manifest['checkpoint_hashes_matched'] == len(ids), 'Export did not check every archive')
    return dict(folder=folder, manifest=manifest, ids=ids, signatures=signatures, matrix=matrix,
                manifest_sha256=digest(folder / 'aggregation_manifest.json'))


def raw_pca(matrix):
    require(matrix.ndim == 2 and min(matrix.shape) >= 2 and np.isfinite(matrix).all(), 'Invalid PCA matrix')
    centered = matrix - matrix.mean(axis=0)
    u, singular, components = np.linalg.svd(centered, full_matrices=False)
    # Same deterministic sign convention as sklearn's full-SVD PCA (svd_flip).
    signs = np.sign(u[np.argmax(np.abs(u), axis=0), np.arange(u.shape[1])])
    signs[signs == 0] = 1
    u *= signs
    components *= signs[:, None]
    variance = singular ** 2 / (len(matrix) - 1)
    require(variance.sum() > 0, 'PCA has no variance')
    return u * singular, components, variance / variance.sum()


def dunn(groups, labels):
    sizes = np.array([len(group) for group in groups])
    values = np.concatenate(groups)
    ranks = stats.rankdata(values)
    _, ties = np.unique(values, return_counts=True)
    n = len(values)
    variance = n * (n + 1) / 12 - np.sum(ties.astype(float) ** 3 - ties) / (12 * (n - 1))
    boundaries = np.r_[0, np.cumsum(sizes)]
    averages = np.array([ranks[boundaries[i]:boundaries[i+1]].mean() for i in range(len(sizes))])
    comparisons = len(groups) * (len(groups) - 1) // 2
    result = []
    for i, j in itertools.combinations(range(len(groups)), 2):
        denominator = np.sqrt(variance * (1 / sizes[i] + 1 / sizes[j]))
        z = (averages[i] - averages[j]) / denominator if denominator > 0 else 0.
        p = float(2 * stats.norm.sf(abs(z)))
        result.append(dict(group1=labels[i], group2=labels[j], n1=int(sizes[i]), n2=int(sizes[j]),
                           z=float(z), p_two_sided=p, p_bonferroni=min(1., comparisons * p)))
    return result


def signature_features(path, signatures):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        channel = reader.fieldnames[0]
        rows = list(reader)
    require(len(rows) == len({r[channel] for r in rows}) == 96, 'Expected 96 unique signature channels')
    result = {}
    for name in signatures:
        p = np.array([float(row[name]) for row in rows])
        require(np.isfinite(p).all() and (p >= 0).all() and p.sum() > 0, 'Invalid signature values')
        p /= p.sum()
        labels = [row[channel] for row in rows]
        require(all(re.fullmatch(r'[ACGT]\[[CT]>[ACGT]\][ACGT]', x) for x in labels), 'Invalid SBS96 label')
        transition = sum(value for label, value in zip(labels, p) if label[2:5] in ('C>T', 'T>C'))
        gc = sum(value for label, value in zip(labels, p) if label[2] == 'C')
        positive = p[p > 0]
        entropy = -float(np.dot(positive, np.log2(positive)))
        result[name] = dict(Signature=name, gc_target_preference=float(gc), entropy_bits=entropy,
                            transition_fraction=float(transition), transversion_fraction=float(1-transition))
    return result


def gc_features(fasta, ids):
    result = {}
    for record in SeqIO.parse(fasta, 'fasta'):
        sequence = str(record.seq[1:-1]).upper()
        require(sequence and len(sequence) % 3 == 0 and set(sequence) <= set('ACGT'), 'Invalid coding sequence')
        identifier = record.id.split('|')[0]
        require(identifier not in result, 'Duplicate FASTA identifier')
        wobble = sequence[2::3]
        result[identifier] = dict(Transcript=identifier, cds_nt=len(sequence),
            gc_total=(sequence.count('G') + sequence.count('C')) / len(sequence),
            gc3=(wobble.count('G') + wobble.count('C')) / len(wobble))
    require(set(result) == set(ids), 'GC feature/transcript table ID mismatch')
    return [result[i] for i in ids]


def correlations(x, y):
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return dict(n=len(x), pearson_r=None, spearman_r=None)
    return dict(n=len(x), pearson_r=float(stats.pearsonr(x, y)[0]), spearman_r=float(stats.spearmanr(x, y)[0]))


def regression(x, y):
    require(len(x) >= 3 and np.ptp(x) > 0, 'Insufficient regression data')
    fit = stats.linregress(x, y)
    return dict(n=len(x), slope=float(fit.slope), intercept=float(fit.intercept),
                r_squared=float(fit.rvalue ** 2), pvalue=float(fit.pvalue), slope_se=float(fit.stderr))


def save_plot(fig, folder, name):
    fig.tight_layout()
    fig.savefig(folder / (name + '.pdf'))
    fig.savefig(folder / (name + '.png'), dpi=180)
    plt.close(fig)


def violin_panel(ax, values, title):
    require(all(len(v) > 1 for v in values), 'Empty/singleton plot group')
    # All observations are used; full tails are retained, with no y-limit censoring.
    for i, array in enumerate(values, 1):
        if np.ptp(array) > 0:
            ax.violinplot([array], positions=[i], showmedians=True, points=100)
        else:
            ax.plot([i-.2, i+.2], [array[0]] * 2)
    ax.set_xticks([1, 2, 3], ['Human', 'Molossus', 'Myotis'], rotation=20)
    ax.set_title(title, fontsize=9)
    ax.set_ylabel('Mean neutral substitutions')


def pairwise(values, label, alternative='two-sided'):
    rows = []
    for i, j in itertools.combinations(range(3), 2):
        test = stats.mannwhitneyu(values[i], values[j], alternative=alternative,
                                 method='asymptotic', use_continuity=True)
        rows.append(dict(comparison=label, group1=SPECIES[i], group2=SPECIES[j], n1=len(values[i]),
            n2=len(values[j]), alternative=alternative, U=float(test.statistic),
            rank_biserial=float(2 * test.statistic / (len(values[i]) * len(values[j])) - 1),
            pvalue=float(test.pvalue), p_bonferroni=min(1., float(test.pvalue) * 3)))
    return rows


def analyze(args):
    started = time.perf_counter()
    human = read_export(args.human_export)
    bats = {species: read_export(args.bat_export / species) for species in SPECIES[1:]}
    enzymes = {species: read_export(args.apobec_export / species) for species in SPECIES}
    datasets = [human, *bats.values(), *enzymes.values()]
    for data in datasets:
        cfg = data['manifest']['configuration']
        require(cfg['criterion'] == 'protein_identity' and cfg['n_trajectories'] == args.n_traj,
                'Unexpected stop criterion or number of paths')
    require(all(enzymes[s]['manifest']['configuration']['mutation_model'] == 'experimental_apobec_NNC_v1'
                for s in SPECIES), 'Expected experimental enzyme runs')
    if not args.allow_subset:
        require(len(human['ids']) == COUNTS[0] and len(human['signatures']) == 79, 'Expected full human dataset')
        require(all(len(enzymes[s]['ids']) == count for s, count in zip(SPECIES, COUNTS)), 'Wrong full enzyme populations')
        require(all(len(bats[s]['ids']) == count and len(bats[s]['signatures']) == 8
                    for s, count in zip(SPECIES[1:], COUNTS[1:])), 'Wrong full bat COSMIC populations')
    require(enzymes['human']['ids'] == human['ids'] and all(enzymes[s]['ids'] == bats[s]['ids'] for s in bats),
            'Transcript IDs/order differ between COSMIC and enzyme analyses')
    require(all(enzymes[s]['manifest']['input_sha256'] == data['manifest']['input_sha256']
                for s, data in dict(human=human, **bats).items()), 'Different DNA inputs across process comparisons')
    require(all(enzymes[s]['signatures'] == enzymes['human']['signatures'] for s in SPECIES), 'Different enzyme columns')
    fasta = Path(human['manifest']['input_path'])
    require(digest(fasta) == human['manifest']['input_sha256'], 'Human FASTA changed since export')
    output = args.out_dir.resolve()
    require(not output.exists(), 'Output directory already exists')
    protected = [data['folder'] for data in datasets] + [fasta.parent.resolve()]
    require(all(output != path and path not in output.parents for path in protected), 'Output nested inside input data')
    signature_sha, fasta_sha = digest(args.signatures), digest(fasta)
    require(all(data['manifest']['configuration']['signature_sha256'] == signature_sha
                for data in [human, *bats.values()]), 'Feature signature file differs from simulation input')
    features = signature_features(args.signatures, human['signatures'])
    gc = gc_features(fasta, human['ids'])
    groups_path = HERE.parent / '02_run_sinabro/manuscript_signature_groups.json'
    groups = dict(all=human['signatures'], **json.loads(groups_path.read_text())['groups'])
    if args.allow_subset:
        groups = {k: [s for s in v if s in human['signatures']] for k, v in groups.items()}
    require(all(len(v) >= 2 for v in groups.values()), 'Insufficient signatures in a manuscript subset')
    report = dict(analysis_passed=False, original_paper_postprocessing=True,
        includes_new_reviewer_analyses=False, bounded_subset=args.allow_subset,
        input_manifests={str(d['folder']): d['manifest_sha256'] for d in datasets},
        signature_sha256=signature_sha, human_fasta_sha256=fasta_sha, groups_sha256=digest(groups_path),
        script_sha256=digest(__file__), numpy=np.__version__, scipy=scipy.__version__, pca={}, statistics={})
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.paper-postprocessing-', dir=output.parent) as temp:
        stage = Path(temp) / 'results'
        stage.mkdir()
        write_csv(stage / 'human_sequence_features.csv', gc)
        write_csv(stage / 'signature_features.csv', features.values())
        human_summary = []
        for group, names in groups.items():
            selected = human['matrix'][:, [human['signatures'].index(name) for name in names]]
            scores, components, ratios = raw_pca(selected)
            write_csv(stage / f'pca_{group}_variance.csv',
                [dict(PC=i+1, explained_variance_ratio=value, cumulative=float(ratios[:i+1].sum())) for i, value in enumerate(ratios)])
            write_csv(stage / f'pca_{group}_loadings.csv', [dict(PC=i+1, Signature=name, coefficient=value)
                for i, row in enumerate(components) for name, value in zip(names, row)])
            write_csv(stage / f'pca_{group}_scores.csv', [dict(Transcript=identifier,
                **{f'PC{i+1}': value for i, value in enumerate(row[:10])}) for identifier, row in zip(human['ids'], scores)])
            report['pca'][group] = dict(signatures=len(names), PC1_PC2_variance=float(ratios[:2].sum()),
                PC1_PC10_variance=float(ratios[:10].sum()),
                pc1_gc3=correlations(scores[:, 0], [r['gc3'] for r in gc]),
                largest_positive_loading={f'PC{i+1}': names[int(np.argmax(row))] for i, row in enumerate(components[:10])})
            samples = [selected[:, i] for i in range(len(names))]
            kw, levene = stats.kruskal(*samples), stats.levene(*samples, center='median')
            report['statistics'][group] = dict(kruskal_Wallis_H=float(kw.statistic), kruskal_p=float(kw.pvalue),
                levene_center='median', levene_F=float(levene.statistic), levene_p=float(levene.pvalue))
            write_csv(stage / f'dunn_bonferroni_{group}.csv', dunn(samples, names))
            means, sds = selected.mean(axis=0), selected.std(axis=0, ddof=1)
            for i, name in enumerate(names):
                human_summary.append(dict(group=group, **features[name], mean_rums=means[i], sd_rums=sds[i]))
            if group == 'noseqerr':
                transition = [features[name]['transition_fraction'] for name in names]
                report['transition_regression'] = regression(transition, means)
                report['entropy_vs_sd'] = correlations([features[name]['entropy_bits'] for name in names], sds)
                fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
                image = axes[0].imshow(components[:10].T, aspect='auto', cmap='coolwarm')
                fig.colorbar(image, ax=axes[0], label='PCA coefficient')
                axes[0].set_yticks(np.arange(len(names)), names, fontsize=4)
                axes[0].set_title('PC1–PC10 loadings')
                axes[1].plot(np.arange(1, len(ratios)+1), ratios.cumsum())
                axes[1].set(xlabel='Principal components', ylabel='Cumulative explained variance')
                scatter = axes[2].scatter(scores[:, 0], scores[:, 1], c=[r['gc3'] for r in gc], s=1, rasterized=True)
                fig.colorbar(scatter, ax=axes[2], label='GC3')
                axes[2].set(xlabel='PC1', ylabel='PC2')
                save_plot(fig, stage, 'figure3_pca')
                fig, axes = plt.subplots(2, 4, figsize=(14, 7))
                for pc, ax in zip(range(2, 10), axes.flat):
                    if pc >= len(components):
                        ax.axis('off')
                        continue
                    order = np.argsort(components[pc])
                    ax.barh(np.arange(len(names)), components[pc][order])
                    ax.set_yticks(np.arange(len(names)), [names[i] for i in order], fontsize=4)
                    ax.set_title(f'PC{pc+1}')
                save_plot(fig, stage, 'figure4_pc_loadings')
                fig, axes = plt.subplots(2, 2, figsize=(9, 7))
                entropy = [features[name]['entropy_bits'] for name in names]
                gc_bias = [features[name]['gc_target_preference'] for name in names]
                for ax, x, y, color, title in [
                    (axes[0,0], entropy, means, gc_bias, 'Mean RUMS vs entropy; GC targeting'),
                    (axes[0,1], entropy, sds, gc_bias, 'SD RUMS vs entropy; GC targeting'),
                    (axes[1,0], entropy, means, transition, 'Mean RUMS vs entropy; transition fraction'),
                    (axes[1,1], transition, means, transition, 'Mean RUMS vs transition fraction')]:
                    ax.scatter(x, y, c=color, s=16)
                    ax.set_title(title, fontsize=9)
                fit = report['transition_regression']
                limits = np.array([min(transition), max(transition)])
                axes[1,1].plot(limits, fit['intercept'] + fit['slope'] * limits)
                save_plot(fig, stage, 'figure5_signature_features')
            print(f'Completed {group}: PCA, features and original rank/variance tests', flush=True)
        write_csv(stage / 'human_signature_feature_summary.csv', human_summary)
        cosmic = dict(human=human, **bats)
        shared = [s for s in ('SBS1','SBS4','SBS6','SBS18','SBS2','SBS13') if all(s in cosmic[d]['signatures'] for d in SPECIES)]
        require(shared, 'No shared COSMIC controls')
        comparisons, distribution_rows = [], []
        fig, axes = plt.subplots(2, 3, figsize=(11, 7))
        for name, ax in zip(shared, axes.flat):
            values = [cosmic[s]['matrix'][:, cosmic[s]['signatures'].index(name)] for s in SPECIES]
            comparisons.extend(pairwise(values, name))
            violin_panel(ax, values, name)
            for species, v in zip(SPECIES, values):
                distribution_rows.append(dict(comparison=name, species=species, n=len(v), mean=float(v.mean()),
                    sd=float(v.std(ddof=1)), median=float(np.median(v)), zero_fraction=float(np.mean(v == 0))))
        for ax in list(axes.flat)[len(shared):]:
            ax.axis('off')
        save_plot(fig, stage, 'figure7_cosmic_controls')
        enzyme_comparisons = []
        fig, axes = plt.subplots(2, 3, figsize=(11, 7))
        for name, ax in zip(enzymes['human']['signatures'], axes.flat):
            values = [enzymes[s]['matrix'][:, enzymes[s]['signatures'].index(name)] for s in SPECIES]
            enzyme_comparisons.extend(pairwise(values, name))
            violin_panel(ax, values, name.replace('HL', ''))
            for species, v in zip(SPECIES, values):
                distribution_rows.append(dict(comparison=name, species=species, n=len(v), mean=float(v.mean()),
                    sd=float(v.std(ddof=1)), median=float(np.median(v)), zero_fraction=float(np.mean(v == 0))))
        save_plot(fig, stage, 'figureS9_experimental_enzymes')
        native = [human['matrix'][:, human['signatures'].index('SBS2')],
                  enzymes['molMol2']['matrix'][:, enzymes['molMol2']['signatures'].index('HLmolMol2_APOBEC3A')],
                  enzymes['myoMyo6']['matrix'][:, enzymes['myoMyo6']['signatures'].index('HLmyoMyo6_APOBEC3A3')]]
        native_rows = pairwise(native, 'native_panel_human_SBS2_bat_APOBEC3A')
        for row, (i,j), direction in zip(native_rows, itertools.combinations(range(3), 2), ('greater','less','less')):
            test = stats.mannwhitneyu(native[i], native[j], alternative=direction,
                                     method='asymptotic', use_continuity=True)
            row.update(original_direction=direction, original_direction_p=float(test.pvalue))
        fig, ax = plt.subplots(figsize=(4, 4))
        violin_panel(ax, native, 'Native panel: human SBS2, bat APOBEC3A')
        save_plot(fig, stage, 'figure7g_native_processes')
        write_csv(stage / 'pairwise_cosmic.csv', comparisons)
        write_csv(stage / 'pairwise_experimental_enzymes.csv', enzyme_comparisons)
        write_csv(stage / 'pairwise_native_processes.csv', native_rows)
        write_csv(stage / 'species_distribution_summary.csv', distribution_rows)
        report['native_panel'] = {s: dict(mean=float(v.mean()), median=float(np.median(v)), sd=float(v.std(ddof=1)))
                                  for s,v in zip(SPECIES,native)}
        report['native_panel_tests'] = native_rows
        require(all(digest(d['folder'] / 'aggregation_manifest.json') == d['manifest_sha256'] for d in datasets)
                and digest(args.signatures) == signature_sha and digest(fasta) == fasta_sha, 'Inputs changed during analysis')
        report.update(analysis_passed=True, elapsed_seconds=time.perf_counter()-started,
            statistical_scope='Original unpaired transcript-level tests; isoforms and repeated signature measures are not independent biological replicates. Not a replacement for revised inferential design.',
            probability_note='Floating-point p=0 indicates numerical underflow, not an exactly zero probability.',
            pca_method='Mean-centered full SVD, no standardization; sklearn svd_flip sign convention; first ten scores exported.',
            remaining_original_analyses=['S1/S2 random-sequence signature validation', 'S11 BLOSUM subset', 'S12 BRCA1 observed-variant overlay'])
        report['outputs'] = [dict(filename=p.name, sha256=digest(p)) for p in sorted(stage.iterdir())]
        (stage / 'analysis_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        require(not output.exists(), 'Output appeared during analysis')
        stage.rename(output)
    print(json.dumps({k:v for k,v in report.items() if k not in ('outputs','input_manifests')},indent=2),flush=True)
    print('Report:', output / 'analysis_summary.json')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--human-export', type=Path, required=True)
    parser.add_argument('--bat-export', type=Path, required=True)
    parser.add_argument('--apobec-export', type=Path, required=True)
    parser.add_argument('--signatures', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--allow-subset', action='store_true', help='Explicit bounded software check, not paper-wide validation')
    analyze(parser.parse_args(argv))


if __name__ == '__main__':
    main()
