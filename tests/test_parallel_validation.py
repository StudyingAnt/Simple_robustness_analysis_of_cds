"""Check that parallel-validation reports reject mismatches and measure overlap."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / '02_run_sinabro'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('parallel_validation_tests_target', SCRIPTS / 'validate_parallel.py')
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.left, self.right = [Path(self.directory.name) / name for name in ('left', 'right')]
        for directory in (self.left, self.right):
            (directory / 'transcripts').mkdir(parents=True)
        self.counts = np.arange(9).reshape(3, 3)
        self.write(self.left)
        self.write(self.right)

    def write(self, directory, counts=None, signatures=None):
        np.savez_compressed(
            directory / 'transcripts/id.npz',
            accepted_steps=self.counts if counts is None else counts,
            terminal_codes=np.zeros((3, 3), dtype=np.uint8),
            metadata=json.dumps(dict(signatures=signatures or ['SBS1', 'SBS2', 'SBS13'],
                                     config_sha256='same-scientific-configuration')))

    def test_equal_per_path_results_are_counted(self):
        self.assertEqual(validation.compare_results(self.left, self.right, 1, 3), 9)

    def test_changed_count_is_rejected(self):
        changed = self.counts.copy()
        changed[0, 0] += 1
        self.write(self.right, counts=changed)
        with self.assertRaises(AssertionError):
            validation.compare_results(self.left, self.right, 1, 3)

    def test_missing_transcript_is_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'Missing'):
            validation.compare_results(self.left, self.right, 2, 3)

    def test_signature_order_mismatch_is_rejected(self):
        self.write(self.right, signatures=['SBS2', 'SBS1', 'SBS13'])
        with self.assertRaisesRegex(AssertionError, 'Signature order'):
            validation.compare_results(self.left, self.right, 1, 3)

    def test_overlapping_case_timestamps_are_measured(self):
        shards = self.left / 'shards'
        shards.mkdir()
        (shards / '000.done.json').write_text(json.dumps(dict(
            case_metrics_file='cases.jsonl', elapsed_seconds=4, compute_and_save_seconds=3,
            workers=2, start_method='spawn', incomplete_paths_in_new_results=0)))
        rows = [dict(started_epoch=0, finished_epoch=2, worker_pid=1),
                dict(started_epoch=1, finished_epoch=3, worker_pid=2)]
        (shards / 'cases.jsonl').write_text('\n'.join(json.dumps(row) for row in rows))
        summary = validation.summarize_execution(self.left)
        self.assertEqual(summary['peak_simultaneous_cases'], 2)
        self.assertEqual(summary['worker_pids'], [1, 2])


if __name__ == '__main__':
    unittest.main()
