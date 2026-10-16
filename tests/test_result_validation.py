"""Reject damaged, incomplete or inconsistent saved results without simulations."""
from contextlib import redirect_stdout
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
import run_rums as runner
import validate_results as validator


class SavedResultTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.make_run()

    def make_run(self, n=3):
        self.root = self.base / f'run_n{n}'
        (self.root / 'transcripts').mkdir(parents=True)
        (self.root / 'shards').mkdir()
        self.fasta = self.base / f'input_n{n}.fa'
        self.fasta.write_text('>id1|gene1\nAATGTAAT\n>id2|gene2\nAATGAAATAAT\n'
                              '>id3|gene3\nAATGAAAAAATAAT\n')
        self.config = dict(runner='rums-v3', input_sha256=validator.digest_file(self.fasta),
                           n_trajectories=n, signatures=['SBS1', 'SBS2'], shards=2,
                           max_steps=10, assignment='greedy_length_balance_v1')
        self.fingerprint = hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()
        (self.root / 'run_manifest.json').write_text(json.dumps(dict(
            configuration=self.config, config_sha256=self.fingerprint, input_path=str(self.fasta))))
        for i in range(1, 4):
            self.save(i)
        for index, assigned in enumerate([1, 2]):
            (self.root / 'shards' / f'{index:03d}.done.json').write_text(json.dumps(dict(
                config_sha256=self.fingerprint, assigned=assigned, skipped=0,
                newly_processed=assigned, incomplete_paths_in_new_results=0)))

    def save(self, index, code=0):
        n = self.config['n_trajectories']
        counts = np.arange(n * 2, dtype=np.int32).reshape(n, 2)
        codes = np.zeros((n, 2), dtype=np.uint8)
        if code:
            codes[0, 1], counts[0, 1] = code, -1
        path = self.root / 'transcripts' / f'id{index}.npz'
        path.unlink(missing_ok=True)
        with patch.multiple(runner, OUTPUT=self.root, N_TRAJ=n, SIGNATURES=['SBS1', 'SBS2'],
                            FINGERPRINT=self.fingerprint, create=True):
            runner.save_transcript(f'id{index}', f'id{index}|gene{index}', dict(
                counts=counts, codes=codes, seconds=np.array([0.1, 0.2])))

    def change(self, **changes):
        path = self.root / 'transcripts/id1.npz'
        with np.load(path, allow_pickle=False) as archive:
            saved = {name: archive[name] for name in archive.files}
        saved.update(changes)
        np.savez_compressed(path, **saved)

    def audit(self):
        with redirect_stdout(io.StringIO()):
            return validator.validate_run(self.root, self.fasta, self.base / 'audit', expected_records=3)

    def test_valid_run_checks_every_path_and_preserves_all_inputs(self):
        before = {str(p): validator.digest_file(p) for p in self.root.rglob('*') if p.is_file()}
        summary = self.audit()
        self.assertTrue(summary['validation_passed'])
        self.assertEqual(summary['checked_paths'], 18)
        self.assertEqual(summary['completed_shards'], 2)
        self.assertEqual(summary['structurally_valid_files'], 3)
        self.assertEqual(before, {str(p): validator.digest_file(p)
                                 for p in self.root.rglob('*') if p.is_file()})

    def test_spawn_cli_checks_results_without_sinabro_or_scheduler(self):
        report = self.base / 'spawn_audit'
        result = subprocess.run([sys.executable, str(SCRIPTS / 'validate_results.py'),
                                 '--run-dir', str(self.root), '--report-dir', str(report),
                                 '--workers', '2'], capture_output=True, text=True, cwd=self.base)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(json.loads((report / 'validation_summary.json').read_text())['validation_passed'])

    def test_incorrect_mean_is_detected(self):
        self.change(mean=np.array([999., 3.]))
        self.assertEqual(self.audit()['invalid_files'], 1)
        self.assertIn('mean does not match', (self.base / 'audit/transcript_audit.csv').read_text())

    def test_signature_order_mismatch_is_detected(self):
        path = self.root / 'transcripts/id1.npz'
        with np.load(path, allow_pickle=False) as saved:
            metadata = json.loads(saved['metadata'].item())
        metadata['signatures'].reverse()
        self.change(metadata=json.dumps(metadata))
        self.assertEqual(self.audit()['invalid_files'], 1)

    def test_missing_and_extra_files_are_detected_by_input_identity(self):
        path = self.root / 'transcripts/id1.npz'
        path.rename(path.with_name('unexpected.npz'))
        summary = self.audit()
        self.assertFalse(summary['validation_passed'])
        self.assertEqual((summary['missing_transcripts'], summary['extra_transcripts']), (1, 1))

    def test_truncated_archive_is_reported_and_other_files_are_checked(self):
        path = self.root / 'transcripts/id1.npz'
        path.write_bytes(path.read_bytes()[:80])
        summary = self.audit()
        self.assertEqual((summary['invalid_files'], summary['structurally_valid_files']), (1, 2))

    def test_unknown_terminal_code_is_rejected(self):
        self.save(1, code=255)
        self.assertEqual(self.audit()['invalid_files'], 1)

    def test_incomplete_results_are_rejected_even_if_resumed_shards_report_zero(self):
        self.save(1, code=2)
        shard = self.root / 'shards/001.done.json'
        data = json.loads(shard.read_text())
        data.update(skipped=2, newly_processed=0)
        shard.write_text(json.dumps(data))
        summary = self.audit()
        self.assertFalse(summary['validation_passed'])
        self.assertEqual(summary['invalid_files'], 0)
        self.assertEqual(summary['max_steps_paths'], 1)
        self.assertFalse(summary['all_paths_complete'])

    def test_single_trajectory_nan_sd_is_valid(self):
        self.make_run(n=1)
        self.assertTrue(self.audit()['validation_passed'])

    def test_wrong_aggregate_complete_count_is_rejected(self):
        self.change(n_complete=np.array([2, 3]))
        self.assertEqual(self.audit()['invalid_files'], 1)

    def test_missing_shard_marker_is_rejected(self):
        (self.root / 'shards/001.done.json').unlink()
        summary = self.audit()
        self.assertFalse(summary['validation_passed'])
        self.assertEqual(summary['completed_shards'], 1)

    def test_mismatched_input_hash_writes_failure_report(self):
        other = self.base / 'other.fa'
        other.write_text(self.fasta.read_text() + '\n')
        summary = validator.validate_run(self.root, other, self.base / 'audit')
        self.assertFalse(summary['validation_passed'])
        self.assertIn('SHA256 differs', summary['issues'][0])
        self.assertTrue((self.base / 'audit/validation_summary.json').is_file())

    def test_report_cannot_overwrite_a_previous_audit_or_enter_run_directory(self):
        self.assertTrue(self.audit()['validation_passed'])
        with self.assertRaises(FileExistsError):
            self.audit()
        with self.assertRaisesRegex(ValueError, 'outside'):
            validator.validate_run(self.root, self.fasta, self.root / 'audit')


if __name__ == '__main__':
    unittest.main()
