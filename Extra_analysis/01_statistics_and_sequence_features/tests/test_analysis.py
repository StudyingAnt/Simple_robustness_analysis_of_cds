"""Independent numerical oracles and end-to-end, read-only workflow checks."""
from contextlib import redirect_stdout
import csv
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest
import warnings

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import statistics_utils as su
import run_analysis as run


class NumericalTests(unittest.TestCase):
    def test_signed_rank_matches_scipy_with_ties_zeros_and_sign_reversal(self):
        for d in (np.array([0, 1, 1, -2, 3, 5, -5, 2, 8, 3, -1, 0]),
                  np.arange(1, 101), np.ones(20), np.array([1, -1]*20)):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                reference = stats.wilcoxon(d, zero_method='wilcox', correction=False, mode='approx')
            result = su.signed_rank(d)
            self.assertAlmostEqual(result['W_min'], reference.statistic)
            self.assertAlmostEqual(result['p'], reference.pvalue, places=13)
            reverse = su.signed_rank(-d)
            self.assertAlmostEqual(result['rank_biserial'], -reverse['rank_biserial'])
            self.assertAlmostEqual(result['p'], reverse['p'])

    def test_all_zeros_and_underflow_are_explicit(self):
        zero = su.signed_rank(np.zeros(30))
        self.assertEqual((zero['p'], zero['rank_biserial'], zero['n_nonzero']), (1, 0, 0))
        tiny = su.signed_rank(np.arange(1, 100001))
        self.assertTrue(tiny['p_underflow'])
        self.assertTrue(np.isfinite(tiny['p_log10']))
        self.assertLess(tiny['p_log10'], -1000)
        json.dumps(tiny, allow_nan=False)

    def test_holm_matches_known_values_and_retains_tiny_order(self):
        p = np.array([.01, .04, .03, .002])
        np.testing.assert_allclose(np.exp(su.holm_log_adjust(np.log(p))), [.03, .06, .06, .008])
        adjusted = su.holm_log_adjust(np.array([-10000., -5000., 0.]))
        self.assertLess(adjusted[0], adjusted[1])
        self.assertEqual(adjusted[2], 0.)

    def test_median_interval_exact_binomial_coverage_and_small_n(self):
        # For n=20, [X_(6), X_(15)] has coverage 1 - 2*P(Bin(20,.5)<=5).
        result = su.median_interval(np.arange(1, 21))
        self.assertEqual((result['ci_low'], result['ci_high']), (6, 15))
        coverage = sum(stats.binom.pmf(b, 20, .5) for b in range(6, 15))
        self.assertAlmostEqual(result['ci_coverage_bound'], coverage)
        self.assertGreaterEqual(coverage, .95)
        self.assertEqual(su.median_interval([1, 2, 3])['ci_low'], None)
        self.assertEqual(su.median_interval(np.zeros(30))['ci_high'], 0)

    def test_paired_median_is_not_difference_of_marginal_medians(self):
        x = np.array([0, 5, 100]*10)
        y = np.array([4, 6, 0]*10)
        matrix = np.column_stack([x, y, x+2])
        _, _, rows = su.paired_analysis(matrix, ['a', 'b', 'c'], 1)
        row = rows[0]
        self.assertEqual(row['paired_median_difference'], -1)
        self.assertEqual(row['median_group1'] - row['median_group2'], 1)

    def test_count_lattice_and_friedman_oracle(self):
        matrix = np.random.default_rng(4).integers(0, 100, (30, 8)) / 1000
        omnibus, _, rows = su.paired_analysis(matrix, list('abcdefgh'), 1000)
        oracle = stats.friedmanchisquare(*np.rint(matrix*1000).T)
        self.assertAlmostEqual(omnibus['statistic'], oracle.statistic)
        self.assertEqual(len(rows), 28)
        self.assertTrue(all(row['holm_family_size'] == 28 for row in rows))
        with self.assertRaisesRegex(ValueError, 'lattice'):
            su.to_count_sums(np.array([[.12345]]), 1000)
        flat, _, _ = su.paired_analysis(np.ones((30, 3)), list('abc'), 1)
        self.assertEqual((flat['statistic'], flat['p']), (0, 1))

    def test_independent_rank_effect_matches_all_pair_signs(self):
        x, y = np.array([1, 3, 3]), np.array([2, 3, 5, 5])
        reference = np.mean([np.sign(a-b) for a, b in itertools.product(x, y)])
        self.assertAlmostEqual(su.independent_rank_effect(x, y), reference)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.names = ['SBS1', 'SBS2', 'SBS13']
        self.enzymes = ['HLmolMol2_APOBEC3A', 'HLmolMol2_APOBEC3C', 'HLmolMol2_APOBEC1',
                        'HLmyoMyo6_APOBEC3A3', 'HLmyoMyo6_APOBEC3C', 'HLmyoMyo6_APOBEC1']
        for species in run.SPECIES:
            fasta = self.root / 'sequences' / species / 'input.fa'
            fasta.parent.mkdir(parents=True)
            fasta.write_text(''.join(f'>{species}{i}\nCATG'+('AAA', 'GGC', 'TTC')[i % 3]*(i+1)+'TAAG\n'
                                     for i in range(30)))
            self.make_export(self.root/'human' if species == 'human' else self.root/'bats'/species,
                             species, fasta, self.names)
            self.make_export(self.root/'enzymes'/species, species, fasta, self.enzymes)

    def make_export(self, folder, species, fasta, names):
        folder.mkdir(parents=True)
        matrix = np.random.default_rng(len(names) + run.SPECIES.index(species)).integers(0, 100, (30, len(names))) / 10
        path = folder / 'rums_mean_by_transcript_id.csv'
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['Transcript', *names])
            writer.writerows([[species+str(i), *row] for i, row in enumerate(matrix)])
        annotations = [dict(transcript_id=species+str(i), cds_nt=3*(i+3), internal_stop='False',
                            **(dict(gene_id=f'gene{i//2}.1', gene_name=f'G{i//2}') if species == 'human' else {}))
                       for i in range(30)]
        run.write_csv(folder / 'transcript_annotations.csv', annotations)
        config = dict(signatures=names, criterion='protein_identity', n_trajectories=10,
                      input_sha256=run.digest(fasta), signature_sha256='test-cosmic', profiles_sha256='test-profile',
                      master_seed=42, max_steps=10000, package_source_sha256='test-package',
                      mutation_model='experimental_apobec_NNC_v1' if names == self.enzymes else 'COSMIC')
        manifest = dict(aggregation_passed=True, configuration=config,
            config_sha256=hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest(),
            transcripts=30, checkpoint_hashes_matched=30, input_path=str(fasta), input_sha256=run.digest(fasta),
            outputs=[dict(filename=p.name, sha256=run.digest(p)) for p in sorted(folder.iterdir())])
        (folder / 'aggregation_manifest.json').write_text(json.dumps(manifest))

    def invoke(self, *extra):
        with redirect_stdout(io.StringIO()):
            return run.main(['--human-export', str(self.root/'human'), '--out-dir', str(self.root/'out'),
                             '--n-traj', '10', '--allow-subset', *extra])

    def test_full_workflow_keeps_inputs_and_separates_species_inference(self):
        before = {p: run.digest(p) for p in self.root.rglob('*') if p.is_file()}
        result = self.invoke('--bat-export', str(self.root/'bats'), '--apobec-export', str(self.root/'enzymes'))
        self.assertTrue(result['analysis_passed'])
        self.assertFalse(result['fresh_simulation'])
        self.assertEqual(result['human_genes'], 15)
        self.assertEqual(result['paired_comparisons'], 3)
        self.assertEqual(result['sequence_feature_datasets'], list(run.SPECIES))
        self.assertEqual(result['species_analysis']['comparisons'], 30)
        self.assertFalse(result['species_analysis']['species_confidence_intervals_completed'])
        self.assertTrue((self.root/'out/gc3_cohort_distributions.pdf').read_bytes().startswith(b'%PDF'))
        with (self.root/'out/human_representative_transcripts.csv').open() as stream:
            selected = list(csv.DictReader(stream))
        self.assertEqual({r['Transcript'] for r in selected}, {'human'+str(i) for i in range(1, 30, 2)})
        self.assertTrue(all(run.digest(p) == sha for p, sha in before.items()))
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.invoke()

    def test_tampering_wrong_population_and_relocated_fasta(self):
        with self.assertRaisesRegex(ValueError, '79'):
            run.main(['--human-export', str(self.root/'human'), '--out-dir', str(self.root/'out'), '--n-traj', '10'])
        data = run.load_dataset(self.root/'human')
        relocated = self.root/'relocated.fa'
        relocated.write_bytes(data['fasta'].read_bytes())
        self.assertEqual(run.load_dataset(self.root/'human', relocated)['ids'], data['ids'])
        path = self.root/'human/transcript_annotations.csv'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'Changed export'):
            self.invoke()
        self.assertFalse((self.root/'out').exists())

    def test_representative_selection_never_depends_on_rums(self):
        data = run.load_dataset(self.root/'human')
        selected = run.representatives(data)
        data['matrix'][:] = np.random.default_rng(99).normal(size=data['matrix'].shape)
        self.assertEqual(run.representatives(data), selected)
        data['annotations'][0]['cds_nt'] = data['annotations'][1]['cds_nt']
        self.assertIn(0, run.representatives(data))
        data['annotations'][0]['gene_id'] = ''
        with self.assertRaisesRegex(ValueError, 'gene_id'):
            run.representatives(data)


if __name__ == '__main__':
    unittest.main()
