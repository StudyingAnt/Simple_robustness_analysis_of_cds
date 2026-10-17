"""Verify descriptive ranks, input correspondence and full-data plot behavior."""
from contextlib import redirect_stdout
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '03_mean_std_analysis'))
import analyze_signature_distributions as analysis
import test_result_aggregation as fixtures


class DescriptiveAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.AggregationTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        groups = self.fixture.base / 'groups.json'
        groups.write_text(json.dumps(dict(groups=dict(noseqerr=['SBS2', 'SBS1'], aetiology=['SBS2', 'SBS1']))))
        self.fixture.export(signature_groups=groups)
        self.tables = self.fixture.output
        self.output = self.fixture.base / 'figures'

    def run_analysis(self):
        with redirect_stdout(io.StringIO()):
            return analysis.analyze(self.tables, self.output, density_points=32, dpi=60, page_size=2)

    def rewrite_inventory_hash(self, name):
        path = self.tables / 'aggregation_manifest.json'
        manifest = json.loads(path.read_text())
        for entry in manifest['outputs']:
            if entry['filename'] == name:
                entry['sha256'] = analysis.digest_file(self.tables / name)
        path.write_text(json.dumps(manifest))

    def test_statistics_distinguish_path_and_transcript_variability_and_ties(self):
        data = np.array([[2., 3., 2.], [1., 0., 1.], [3., 2., 3.]])
        rows = analysis.signature_statistics(['SBS2', 'SBS1', 'SBS3'], data)
        self.assertEqual(rows[0]['mean_rums_across_transcripts'], 2.)
        self.assertEqual(rows[0]['sd_rums_across_transcripts'], 1.)
        ranks = analysis.rank_groups(rows, dict(all=['SBS2', 'SBS1', 'SBS3']))
        self.assertEqual([r['Signature'] for r in ranks if r['ranking_metric'] == 'mean'],
                         ['SBS2', 'SBS3', 'SBS1'])
        self.assertEqual([r['Signature'] for r in ranks if r['ranking_metric'] == 'sd'],
                         ['SBS1', 'SBS2', 'SBS3'])

    def test_density_includes_full_support_and_handles_constant_columns(self):
        _, constant = analysis.density((0, np.zeros(7), 32))
        self.assertEqual(constant['grid'].tolist(), [0.])
        values = np.array([0., 0., 1., 2., 100.])
        _, curve = analysis.density((0, values, 32))
        self.assertEqual((curve['grid'][0], curve['grid'][-1]), (0., 100.))
        self.assertAlmostEqual(curve['relative_density'].max(), 1.)

    def test_modified_csv_is_rejected_before_plotting(self):
        path = self.tables / 'rums_mean_by_transcript_id.csv'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'Changed aggregation'):
            self.run_analysis()
        self.assertFalse(self.output.exists())

    def test_incorrect_rehashed_summary_is_detected(self):
        name = 'signature_summary.csv'
        path = self.tables / name
        with path.open(newline='') as stream:
            reader = csv.DictReader(stream)
            header, rows = reader.fieldnames, list(reader)
        rows[0]['mean_rums_across_transcripts'] = '999'
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=header)
            writer.writeheader()
            writer.writerows(rows)
        self.rewrite_inventory_hash(name)
        with self.assertRaisesRegex(ValueError, 'summary mismatch'):
            self.run_analysis()

    def test_reordered_legacy_rows_are_rejected_even_with_valid_hash(self):
        name = 'all_gencode_rums_profile.csv'
        path = self.tables / name
        lines = path.read_text().splitlines()
        path.write_text('\n'.join([lines[0], *reversed(lines[1:])]) + '\n')
        self.rewrite_inventory_hash(name)
        with self.assertRaisesRegex(ValueError, 'row correspondence'):
            self.run_analysis()

    def test_spawn_cli_renders_pdf_png_and_records_observed_leaders(self):
        before = {p.name: analysis.digest_file(p) for p in self.tables.iterdir()}
        result = subprocess.run([sys.executable, str(ROOT / '03_mean_std_analysis/analyze_signature_distributions.py'),
                                 '--aggregation-dir', str(self.tables), '--out-dir', str(self.output),
                                 '--workers', '2', '--density-points', '32', '--dpi', '60', '--page-size', '2'],
                                capture_output=True, text=True, cwd=self.fixture.base)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = json.loads((self.output / 'analysis_summary.json').read_text())
        self.assertTrue(summary['analysis_passed'])
        self.assertTrue(summary['summary_and_rank_agreement'])
        self.assertFalse(summary['subsampling'])
        self.assertEqual(summary['leaders']['aetiology']['mean']['signature'], 'SBS2')
        self.assertEqual(summary['leaders']['aetiology']['sd']['signature'], 'SBS1')
        self.assertEqual(summary['manuscript_claim_checks']['aetiology']['SBS2_sd_rank'], 2)
        self.assertEqual(summary['plotting']['pdf_rank_pages'], 6)
        self.assertEqual((self.output / 'aetiology_top_bottom.pdf').read_bytes()[:4], b'%PDF')
        self.assertEqual((self.output / 'aetiology_top_bottom.png').read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
        self.assertEqual(before, {p.name: analysis.digest_file(p) for p in self.tables.iterdir()})
        for entry in summary['outputs']:
            self.assertEqual(entry['sha256'], analysis.digest_file(self.output / entry['filename']))

    def test_existing_or_nested_output_and_wrong_expected_dimensions_are_rejected(self):
        self.output.mkdir()
        with self.assertRaisesRegex(ValueError, 'exists'):
            self.run_analysis()
        with self.assertRaisesRegex(ValueError, 'outside'):
            analysis.analyze(self.tables, self.tables / 'nested')
        with self.assertRaisesRegex(ValueError, 'transcript count'):
            analysis.analyze(self.tables, self.fixture.base / 'other', expected_records=4)


if __name__ == '__main__':
    unittest.main()
