"""Check exported identifiers, audit binding and distinct levels of variation."""
import ast
from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / '02_run_sinabro'
sys.path.insert(0, str(SCRIPTS))
import aggregate_results as exporter
import run_rums as runner
import validate_results as validator


class AggregationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root, self.audit, self.output = [self.base / name for name in ('run', 'audit', 'tables')]
        (self.root / 'transcripts').mkdir(parents=True)
        (self.root / 'shards').mkdir()
        self.fasta = self.base / 'input.fa'
        headers = [f'ENST_test{i}|ENSG_test{1 if i < 3 else 2}|OTTHG|OTTHT|'
                   f'GENE{1 if i < 3 else 2}-{201 if i != 2 else 202}|GENE{1 if i < 3 else 2}|'
                   for i in range(1, 4)]
        sequences = ['AATGAAATAAT', 'AATGTAATAAT', 'AATGTAAT']
        self.fasta.write_text(''.join(f'>{headers[i]}\n{sequences[i]}\n' for i in (1, 0, 2)))
        self.config = dict(runner='rums-v3', input_sha256=validator.digest_file(self.fasta),
                           n_trajectories=3, signatures=['SBS2', 'SBS1'], shards=1,
                           max_steps=10, assignment='greedy_length_balance_v1')
        self.fingerprint = hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()
        (self.root / 'run_manifest.json').write_text(json.dumps(dict(
            config_sha256=self.fingerprint, configuration=self.config, input_path=str(self.fasta))))
        with patch.multiple(runner, OUTPUT=self.root, N_TRAJ=3, SIGNATURES=['SBS2', 'SBS1'],
                            FINGERPRINT=self.fingerprint, create=True):
            for i, counts in enumerate([[[0, 1], [2, 3], [4, 5]],
                                       [[1, 0]] * 3, [[3, 2]] * 3]):
                runner.save_transcript(f'ENST_test{i+1}', headers[i], dict(
                    counts=np.array(counts, dtype=np.int32), codes=np.zeros((3, 2), dtype=np.uint8),
                    seconds=np.array([0.1, 0.2])))
        (self.root / 'shards/000.done.json').write_text(json.dumps(dict(
            config_sha256=self.fingerprint, assigned=3, skipped=0, newly_processed=3,
            incomplete_paths_in_new_results=0)))
        with redirect_stdout(io.StringIO()):
            self.assertTrue(validator.validate_run(self.root, self.fasta, self.audit)['validation_passed'])

    def export(self, **kwargs):
        with redirect_stdout(io.StringIO()):
            return exporter.aggregate(self.root, self.audit, self.fasta, self.output, **kwargs)

    def rows(self, name):
        with (self.output / name).open(newline='') as stream:
            return list(csv.DictReader(stream))

    def test_export_preserves_values_identity_and_isoforms(self):
        before = {str(p): validator.digest_file(p) for folder in (self.root, self.audit)
                  for p in folder.rglob('*') if p.is_file()}
        result = self.export()
        self.assertTrue(result['aggregation_passed'])
        self.assertEqual((result['transcripts'], result['unique_genes']), (3, 2))
        self.assertEqual(result['internal_stop_transcripts_retained'], 1)
        self.assertEqual(result['total_trajectories'], 18)
        canonical = self.rows('rums_mean_by_transcript_id.csv')
        legacy = self.rows('all_gencode_rums_profile.csv')
        self.assertEqual(list(canonical[0]), ['Transcript', 'SBS2', 'SBS1'])
        self.assertEqual([r['Transcript'] for r in canonical], ['ENST_test1', 'ENST_test2', 'ENST_test3'])
        self.assertEqual([r['Transcript'] for r in legacy], ['GENE1-201', 'GENE1-202', 'GENE2-201'])
        self.assertEqual([float(r['SBS2']) for r in canonical], [2., 1., 3.])
        self.assertEqual([float(r['SBS1']) for r in canonical], [3., 0., 2.])
        self.assertEqual([r['SBS2'] for r in canonical], [r['SBS2'] for r in legacy])
        self.assertEqual(float(self.rows('rums_sd_by_transcript_id.csv')[0]['SBS2']), 2.)
        self.assertAlmostEqual(float(self.rows('rums_mc_se_by_transcript_id.csv')[0]['SBS2']), 2/np.sqrt(3))
        summaries = self.rows('signature_summary.csv')
        self.assertEqual(float(summaries[0]['sd_rums_across_transcripts']), 1.)
        self.assertEqual(summaries[0]['rank_mean_desc'], '1')
        self.assertEqual(summaries[1]['rank_sd_desc'], '1')
        self.assertEqual(before, {str(p): validator.digest_file(p) for folder in (self.root, self.audit)
                                 for p in folder.rglob('*') if p.is_file()})
        for entry in result['outputs']:
            self.assertEqual(entry['sha256'], validator.digest_file(self.output / entry['filename']))

    def test_optional_subset_order_is_preserved(self):
        groups = self.base / 'groups.json'
        groups.write_text(json.dumps(dict(groups=dict(noseqerr=['SBS1'], aetiology=['SBS2']))))
        result = self.export(signature_groups=groups)
        self.assertEqual(result['signature_groups'], dict(noseqerr=['SBS1'], aetiology=['SBS2']))
        row = self.rows('all_gencode_rums_profile_noseqerr.csv')[0]
        self.assertEqual(row, dict(Transcript='GENE1-201', SBS1='3.0'))

    def test_changed_checkpoint_is_rejected_without_publishing_partial_tables(self):
        path = self.root / 'transcripts/ENST_test2.npz'
        path.write_bytes(path.read_bytes() + b'changed after audit')
        with self.assertRaisesRegex(ValueError, 'changed since validation'):
            self.export()
        self.assertFalse(self.output.exists())

    def test_failed_audit_cannot_be_exported(self):
        path = self.audit / 'validation_summary.json'
        value = json.loads(path.read_text())
        value['validation_passed'] = False
        path.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, 'passing'):
            self.export()

    def test_mismatched_manifest_is_rejected(self):
        path = self.root / 'run_manifest.json'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'different simulation'):
            self.export()

    def test_missing_checkpoint_and_inventory_entry_are_rejected(self):
        path = self.audit / 'transcript_audit.csv'
        lines = path.read_text().splitlines()
        path.write_text('\n'.join(lines[:-1]) + '\n')
        with self.assertRaisesRegex(ValueError, 'cover the input exactly'):
            self.export()
        (self.root / 'transcripts/ENST_test1.npz').unlink()
        with self.assertRaisesRegex(ValueError, 'file set changed'):
            self.export()

    def test_input_expected_count_and_unavailable_subsets_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'transcript count'):
            self.export(expected_records=4)
        groups = self.base / 'groups.json'
        groups.write_text(json.dumps(dict(groups=dict(noseqerr=['SBS99']))))
        with self.assertRaisesRegex(ValueError, 'Unavailable'):
            self.export(signature_groups=groups)

    def test_existing_output_and_protected_directories_are_rejected(self):
        self.export()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.export()
        with self.assertRaisesRegex(ValueError, 'outside'):
            exporter.aggregate(self.root, self.audit, self.fasta, self.root / 'tables')

    def test_spawn_cli_works_without_sinabro_or_scheduler(self):
        result = subprocess.run([sys.executable, str(SCRIPTS / 'aggregate_results.py'),
                                 '--run-dir', str(self.root), '--validation-dir', str(self.audit),
                                 '--out-dir', str(self.output), '--workers', '2'],
                                capture_output=True, text=True, cwd=self.base)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(json.loads((self.output / 'aggregation_manifest.json').read_text())['aggregation_passed'])


class ManuscriptGroupsTests(unittest.TestCase):
    def test_group_membership_and_order_match_preserved_manuscript_code(self):
        data = json.loads((SCRIPTS / 'manuscript_signature_groups.json').read_text())
        original = SCRIPTS.parent / data['source']
        self.assertEqual(data['source_sha256'], validator.digest_file(original))
        for script in (original, SCRIPTS.parent / 'historical/04_PCA_and_gc_contents/run_pca.py'):
            for node in ast.parse(script.read_text()).body:
                if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                    name = node.targets[0].id
                    if name in ('signatures_noseqerr', 'signatures_aetiology'):
                        values = ast.literal_eval(node.value)
                        if values[0] == 'Transcript':
                            values = values[1:]
                        self.assertEqual(values, data['groups'][name.replace('signatures_', '')])
        self.assertEqual(len(data['groups']['noseqerr']), 60)
        self.assertEqual(len(data['groups']['aetiology']), 41)


if __name__ == '__main__':
    unittest.main()
