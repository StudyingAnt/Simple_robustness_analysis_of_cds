"""Additional paired statistics and sequence features from audited RUMS exports.

No mutations are simulated. Existing manuscript workflows and inputs are read-only.
"""
import argparse
import csv
import importlib.util
import itertools
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import numpy as np
import scipy
from scipy import stats
import statistics_utils as paired

# Reuse the existing checked reader and sequence-coordinate definition unchanged.
_spec = importlib.util.spec_from_file_location('sinabro_original_postprocessing',
    REPO / '07_manuscript_validation/validate_manuscript_results.py')
paper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(paper)
require, digest, write_csv = paper.require, paper.digest, paper.write_csv
SPECIES = ('human', 'molMol2', 'myoMyo6')
EXPECTED_COUNTS = dict(zip(SPECIES, (73214, 31185, 23213)))


def load_dataset(folder, fasta_override=None):
    data = paper.read_export(folder)
    cfg = data['manifest']['configuration']
    require(cfg['criterion'] == 'protein_identity', 'Only protein-identity exports are supported')
    required = {'rums_mean_by_transcript_id.csv', 'transcript_annotations.csv'}
    inventory = {row['filename'] for row in data['manifest']['outputs']}
    require(required <= inventory, 'Export must authenticate both mean matrix and annotations')
    with (data['folder'] / 'transcript_annotations.csv').open(newline='') as stream:
        annotations = list(csv.DictReader(stream))
    ids = [row['transcript_id'] for row in annotations]
    require(len(ids) == len(set(ids)) and set(ids) == set(data['ids']), 'Annotation ID mismatch')
    lookup = dict(zip(ids, annotations))
    data['annotations'] = [lookup[name] for name in data['ids']]
    fasta = Path(fasta_override or data['manifest']['input_path']).resolve()
    require(digest(fasta) == data['manifest']['input_sha256'] == cfg['input_sha256'],
            'FASTA differs from audited simulation input')
    data['fasta'] = fasta
    data['features'] = paper.gc_features(fasta, data['ids'])
    for row, feature in zip(data['annotations'], data['features']):
        require(int(row['cds_nt']) == feature['cds_nt'], 'Annotated CDS length disagrees with FASTA')
    paired.to_count_sums(data['matrix'], cfg['n_trajectories'])
    return data


def representatives(data):
    selected = {}
    for index, row in enumerate(data['annotations']):
        require(row.get('gene_id', '').strip(), 'Human representative selection requires gene_id')
        gene = row['gene_id'].split('.')[0]
        # RUMS values never enter selection. Break equal-length ties by versioned ID.
        key = (-int(row['cds_nt']), row['transcript_id'])
        if gene not in selected or key < selected[gene][0]:
            selected[gene] = (key, index)
    return [selected[gene][1] for gene in sorted(selected)]


def annotated_features(data, species):
    rows = []
    for annotation, feature in zip(data['annotations'], data['features']):
        internal = annotation.get('internal_stop', '')
        if internal == '':
            internal = str(int(annotation.get('internal_stop_count', '0')) > 0)
        rows.append(dict(dataset=species, **feature, gene_id=annotation.get('gene_id', ''),
                         gene_name=annotation.get('gene_name', ''), internal_stop=internal))
    return rows


def correlation(x, y):
    return (None if np.ptp(x) == 0 or np.ptp(y) == 0
            else float(stats.spearmanr(x, y)[0]))


def plot_gc3(features, output):
    fig, ax = paper.plt.subplots(figsize=(6, 4))
    for species, rows in features.items():
        x = np.sort([row['gc3'] for row in rows])
        ax.step(x, np.arange(1, len(x)+1) / len(x), where='post', label=f'{species} (n={len(x):,})')
    ax.set(xlabel='GC fraction at third CDS positions (GC3)', ylabel='Cumulative fraction of transcripts',
           xlim=(0, 1), ylim=(0, 1))
    ax.legend(frameon=False)
    paper.save_plot(fig, output, 'gc3_cohort_distributions')


def species_comparisons(datasets, enzymes, stage):
    summary_rows, contrasts = [], []
    shared = [s for s in datasets['human']['signatures']
              if all(s in datasets[species]['signatures'] for species in SPECIES)]
    require(shared, 'No shared COSMIC signature across cohorts')
    cases = []
    for name in shared:
        cases.append(('fixed_cosmic', name, {s: (datasets[s], name) for s in SPECIES}))
    for name in enzymes['human']['signatures']:
        cases.append(('fixed_experimental', name, {s: (enzymes[s], name) for s in SPECIES}))
    native = dict(human='SBS2', molMol2='HLmolMol2_APOBEC3A', myoMyo6='HLmyoMyo6_APOBEC3A3')
    require(all(name in (datasets['human'] if s == 'human' else enzymes[s])['signatures']
                for s, name in native.items()), 'Missing native-panel profile')
    cases.append(('different_native_profiles', 'native_diagonal',
                  {s: (datasets[s] if s == 'human' else enzymes[s], name) for s, name in native.items()}))
    for kind, label, mapping in cases:
        values = {s: data['matrix'][:, data['signatures'].index(name)]
                  for s, (data, name) in mapping.items()}
        for species, array in values.items():
            summary_rows.append(dict(comparison_type=kind, comparison=label, dataset=species,
                Signature=mapping[species][1], **paired.distribution(array)))
        for left, right in itertools.combinations(SPECIES, 2):
            contrasts.append(dict(comparison_type=kind, comparison=label, group1=left, group2=right,
                profile1=mapping[left][1], profile2=mapping[right][1],
                n1=len(values[left]), n2=len(values[right]),
                difference_of_cohort_medians=float(np.median(values[left]) - np.median(values[right])),
                independent_rank_biserial=paired.independent_rank_effect(values[left], values[right]),
                inferential_pvalue='', ci_low='', ci_high='',
                inference_scope='descriptive annotated transcript cohorts; gene grouping unverified'))
    write_csv(stage / 'species_profile_distributions.csv', summary_rows)
    write_csv(stage / 'species_descriptive_contrasts.csv', contrasts)
    return dict(shared_cosmic_profiles=shared, fixed_experimental_profiles=enzymes['human']['signatures'],
                comparisons=len(contrasts), biological_species_inference=False,
                species_confidence_intervals_completed=False,
                reason='No verified bat gene grouping or independent species replication; no transcript bootstrap CI claimed.')


def analyze(args):
    started = time.perf_counter()
    require(bool(args.bat_export) == bool(args.apobec_export), 'Supply both bat and APOBEC exports, or neither')
    require(args.top_outliers > 0, 'Outlier count must be positive')
    human = load_dataset(args.human_export, args.human_fasta)
    datasets = {'human': human}
    enzymes = {}
    if args.bat_export:
        for species, override in [('molMol2', args.molmol2_fasta), ('myoMyo6', args.myomyo6_fasta)]:
            datasets[species] = load_dataset(args.bat_export / species, override)
        for species in SPECIES:
            data = datasets[species]
            enzymes[species] = load_dataset(args.apobec_export / species, data['fasta'])
            require(enzymes[species]['ids'] == data['ids'], 'COSMIC/enzyme transcript order mismatch')
            require(enzymes[species]['manifest']['configuration']['mutation_model'] == 'experimental_apobec_NNC_v1',
                    'Expected experimental NNC model')
        require(all(enzymes[s]['signatures'] == enzymes['human']['signatures'] for s in SPECIES),
                'Experimental profile order mismatch')
    all_data = list(datasets.values()) + list(enzymes.values())
    require(all(d['manifest']['configuration']['n_trajectories'] == args.n_traj for d in all_data),
            'Unexpected number of trajectories')
    for key in ('package_source_sha256', 'master_seed', 'max_steps'):
        reference = human['manifest']['configuration'].get(key)
        require(reference is not None and all(d['manifest']['configuration'].get(key) == reference for d in all_data),
                'Different or undocumented simulation setting: ' + key)
    require(all(d['manifest']['configuration']['signature_sha256'] == human['manifest']['configuration']['signature_sha256']
                for d in datasets.values()), 'Different COSMIC source tables')
    if enzymes:
        require(all(d['manifest']['configuration']['profiles_sha256'] == enzymes['human']['manifest']['configuration']['profiles_sha256']
                    for d in enzymes.values()), 'Different experimental profile sources')
    if not args.allow_subset:
        require(len(human['signatures']) == 79, 'Full analysis requires 79 human COSMIC profiles')
        for species, data in datasets.items():
            require(len(data['ids']) == EXPECTED_COUNTS[species], 'Unexpected production cohort size')
            if species != 'human':
                require(len(data['signatures']) == 8, 'Full bat COSMIC export requires eight profiles')
        require(not enzymes or len(enzymes['human']['signatures']) == 6, 'Expected six experimental profiles')
    selected = representatives(human)
    output = args.out_dir.resolve()
    protected = [d['folder'] for d in all_data] + [d['fasta'].parent for d in all_data]
    require(not output.exists(), 'Output already exists; choose a new output directory')
    require(all(p != output and p not in output.parents and output not in p.parents for p in protected),
            'Output must be separate from input directories')
    report = dict(analysis_passed=False, stage='additional_statistics_and_sequence_features',
        fresh_simulation=False, bounded_subset=args.allow_subset,
        n_trajectories=args.n_traj, human_transcripts=len(human['ids']), human_genes=len(selected),
        top_outliers_per_signature=args.top_outliers,
        human_signature_count=len(human['signatures']),
        representative_rule='longest retained CDS per stable gene ID; versioned transcript ID lexical tie-break',
        analysis_unit='one representative human CDS per gene; all-transcript distributions descriptive',
        paired_test='two-sided Wilcoxon normal approximation; zeros dropped; average ties; no continuity correction',
        multiple_testing='Holm across all unordered pairs among the supplied human signatures',
        confidence_interval='pointwise 95% binomial/order-statistic CI for median of within-gene differences; not simultaneous',
        signed_rank_caution='Signed-rank inference assumes a symmetric distribution of differences under its null; not a general median test.',
        uncertainty_scope='Conditional sequence/gene sampling summaries; no population/phylogenetic or Monte Carlo-error claim.',
        input_manifests={str(d['folder']): d['manifest_sha256'] for d in all_data},
        fasta_inputs={s: dict(path=str(d['fasta']), sha256=digest(d['fasta'])) for s, d in datasets.items()},
        source_sha256={p.name: digest(p) for p in [Path(__file__), HERE / 'statistics_utils.py',
                                                 REPO / '07_manuscript_validation/validate_manuscript_results.py']},
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        reviewer_coverage={'R2.2.4': 'paired statistics/effects/CI on human genes',
                          'R2.3': 'GC3 and CDS-length distributions; does not standardize orthologues',
                          'R2.6.4': 'annotated robustness outliers',
                          'R2.2.2': 'direct existing-result SBS2/SBS13 comparison',
                          'R1': 'effect/uncertainty/composition reporting; no null model or simulator benchmark claimed'})
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.extra-analysis-', dir=output.parent) as temp:
        stage = Path(temp) / 'output'
        stage.mkdir()
        selected_rows = [dict(gene_id=human['annotations'][i]['gene_id'], Transcript=human['ids'][i],
                              cds_nt=human['features'][i]['cds_nt']) for i in selected]
        write_csv(stage / 'human_representative_transcripts.csv', selected_rows)
        print(f"Paired statistics: {len(selected)} genes x {len(human['signatures'])} signatures", flush=True)
        omnibus, summaries, pairs = paired.paired_analysis(human['matrix'][selected], human['signatures'], args.n_traj,
            progress=lambda n, total: print(f'Paired contrasts {n}/{total}', flush=True))
        report['friedman'] = omnibus
        report['paired_comparisons'] = len(pairs)
        write_csv(stage / 'human_gene_signature_summary.csv', summaries)
        write_csv(stage / 'human_gene_paired_comparisons.csv', pairs)
        transcript_summary = [dict(Signature=s, **paired.distribution(human['matrix'][:, j]))
                              for j, s in enumerate(human['signatures'])]
        write_csv(stage / 'human_transcript_signature_summary.csv', transcript_summary)
        report['sbs2_vs_sbs13'] = next((row for row in pairs if {row['group1'], row['group2']} == {'SBS2', 'SBS13'}), None)
        report['highest_gene_median_signature'] = max(summaries, key=lambda r: r['median'])['Signature']
        features = {s: annotated_features(d, s) for s, d in datasets.items()}
        for species, rows in features.items():
            write_csv(stage / (species + '_sequence_features.csv'), rows)
        feature_summary = [dict(dataset=s, feature=key, **paired.distribution([row[key] for row in rows]))
                           for s, rows in features.items() for key in ('gc3', 'gc_total', 'cds_nt')]
        write_csv(stage / 'sequence_feature_summary.csv', feature_summary)
        plot_gc3(features, stage)
        outliers = []
        for j, name in enumerate(human['signatures']):
            # Export a fixed top N per profile, with deterministic accession tie breaking.
            order = sorted(range(len(human['ids'])), key=lambda i: (-human['matrix'][i, j], human['ids'][i]))
            for rank, i in enumerate(order[:args.top_outliers], 1):
                outliers.append(dict(Signature=name, rank=rank, mean_rums=float(human['matrix'][i, j]),
                                     **features['human'][i]))
        write_csv(stage / 'human_top_outliers_by_signature.csv', outliers)
        correlation_rows = []
        for scope, indices in [('all_transcripts_descriptive', list(range(len(human['ids'])))),
                               ('representative_genes_descriptive', selected)]:
            for j, name in enumerate(human['signatures']):
                for key in ('gc3', 'cds_nt'):
                    x = np.array([features['human'][i][key] for i in indices])
                    correlation_rows.append(dict(scope=scope, Signature=name, feature=key, n=len(indices),
                        spearman_rho=correlation(x, human['matrix'][indices, j])))
        write_csv(stage / 'human_composition_correlations.csv', correlation_rows)
        if enzymes:
            report['species_analysis'] = species_comparisons(datasets, enzymes, stage)
        else:
            report['species_analysis'] = dict(included=False, reason='Human-only invocation')
        report['sequence_feature_datasets'] = list(features)
        # Re-check manifests, all authenticated export files, and FASTAs before publishing.
        for d in all_data:
            reread = paper.read_export(d['folder'])
            require(reread['manifest_sha256'] == d['manifest_sha256'], 'Input manifest changed during analysis')
            require(digest(d['fasta']) == d['manifest']['input_sha256'], 'FASTA changed during analysis')
        report.update(analysis_passed=True, elapsed_seconds=time.perf_counter()-started,
            outputs=[dict(filename=p.name, sha256=digest(p), bytes=p.stat().st_size)
                     for p in sorted(stage.iterdir()) if p.is_file()])
        (stage / 'additional_analysis_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        require(not output.exists(), 'Output appeared during analysis')
        stage.rename(output)
    print(json.dumps({k: report[k] for k in ('analysis_passed', 'fresh_simulation', 'bounded_subset',
        'human_transcripts', 'human_genes', 'human_signature_count', 'paired_comparisons', 'friedman',
        'highest_gene_median_signature', 'sequence_feature_datasets', 'elapsed_seconds')}, indent=2), flush=True)
    print('Report:', output / 'additional_analysis_summary.json', flush=True)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--human-export', type=Path, required=True)
    parser.add_argument('--bat-export', type=Path)
    parser.add_argument('--apobec-export', type=Path)
    parser.add_argument('--human-fasta', type=Path, help='Relocated FASTA; must match recorded hash')
    parser.add_argument('--molmol2-fasta', type=Path)
    parser.add_argument('--myomyo6-fasta', type=Path)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--n-traj', type=int, default=1000)
    parser.add_argument('--top-outliers', type=int, default=10)
    parser.add_argument('--allow-subset', action='store_true', help='Explicit bounded validation, not full-paper inference')
    return analyze(parser.parse_args(argv))


if __name__ == '__main__':
    main()
