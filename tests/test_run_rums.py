"""Check resume safety and CLI validation without generating trajectories."""
from contextlib import redirect_stderr
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / '02_run_sinabro/run_rums.py'
spec = importlib.util.spec_from_file_location('rums_runner', SCRIPT)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        output = Path(self.directory.name)
        (output / 'transcripts').mkdir()
        settings = patch.multiple(
            runner, OUTPUT=output, N_TRAJ=3, SIGNATURES=['SBS1', 'SBS2'],
            FINGERPRINT='checked-configuration', create=True)
        settings.start()
        self.addCleanup(settings.stop)
        self.identifier = 'ENST_test.1'
        self.item = {
            'counts': np.array([[2, 0], [4, -1], [6, 2]], dtype=np.int32),
            'codes': np.array([[0, 0], [0, 2], [0, 0]], dtype=np.uint8),
            'seconds': np.array([1.0, 2.0]),
        }

    def save(self):
        return runner.save_transcript(self.identifier, self.identifier + '|GENE', self.item)

    def check(self, identifier=None):
        return runner.check_checkpoint(runner.checkpoint_path(self.identifier),
                                       identifier or self.identifier)

    def test_missing_checkpoint_is_not_skipped(self):
        self.assertFalse(self.check())

    def test_round_trip_resume_and_incomplete_statistics(self):
        self.assertEqual(self.save(), 1)
        self.assertTrue(self.check())
        with np.load(runner.checkpoint_path(self.identifier), allow_pickle=False) as saved:
            np.testing.assert_array_equal(saved['accepted_steps'], self.item['counts'])
            np.testing.assert_array_equal(saved['terminal_codes'], self.item['codes'])
            np.testing.assert_array_equal(saved['n_complete'], [3, 2])
            self.assertEqual(saved['mean'][0], 4)
            self.assertEqual(saved['median'][0], 4)
            self.assertEqual(saved['sd'][0], 2)
            self.assertAlmostEqual(saved['mc_se'][0], 2 / np.sqrt(3))
            for name in ('mean', 'median', 'sd', 'mc_se'):
                self.assertTrue(np.isnan(saved[name][1]))

    def test_existing_checkpoint_cannot_be_overwritten(self):
        self.save()
        path = runner.checkpoint_path(self.identifier)
        before = path.read_bytes()
        with self.assertRaises(FileExistsError):
            self.save()
        self.assertEqual(path.read_bytes(), before)

    def test_configuration_mismatch_is_rejected(self):
        self.save()
        with patch.object(runner, 'FINGERPRINT', 'different-configuration'):
            with self.assertRaisesRegex(RuntimeError, 'configuration mismatch'):
                self.check()

    def test_wrong_transcript_is_rejected(self):
        self.save()
        with self.assertRaisesRegex(RuntimeError, 'configuration mismatch'):
            self.check('other_transcript')

    def test_unfinished_checkpoint_is_rejected(self):
        self.item['codes'][0, 0] = 255
        self.save()
        with self.assertRaisesRegex(RuntimeError, 'Unfinished checkpoint'):
            self.check()

    def test_invalid_complete_count_is_rejected(self):
        self.item['counts'][0, 0] = -1
        self.save()
        with self.assertRaisesRegex(RuntimeError, 'Invalid checkpoint counts'):
            self.check()


class ArgumentTests(unittest.TestCase):
    def arguments(self):
        return ['--input', 'input.fa', '--out-dir', 'output',
                '--sinabro-source', 'source', '--workers', '8',
                '--n-shards', '100', '--shard-index', '17']

    def test_named_arguments_preserve_manuscript_defaults(self):
        args = runner.parse_arguments(self.arguments())
        self.assertEqual((args.n_traj, args.seed, args.max_steps), (1000, 42, 10000))
        self.assertEqual((args.workers, args.shard_index), (8, 17))

    def test_invalid_parameters_are_rejected_before_execution(self):
        for option, value in [('--workers', '0'), ('--shard-index', '-1'),
                              ('--shard-index', '100'), ('--n-traj', '0'),
                              ('--max-steps', '0'), ('--seed', '-1'), ('--n-shards', '0')]:
            with self.subTest(option=option, value=value), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    runner.parse_arguments(self.arguments() + [option, value])
                self.assertEqual(error.exception.code, 2)

    def test_workstation_defaults_do_not_require_a_scheduler(self):
        args = runner.parse_arguments(['--input', 'input.fa', '--out-dir', 'output'])
        self.assertEqual((args.workers, args.n_shards, args.shard_index), (1, 1, 0))
        self.assertEqual(args.start_method, 'spawn')
        self.assertIsNone(args.sinabro_source)

    def test_duplicate_signatures_are_rejected(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            runner.parse_arguments(self.arguments() + ['--signatures', 'SBS1', 'SBS1'])


class AssignmentAndLockTests(unittest.TestCase):
    def test_shards_cover_input_once_and_preserve_assignment(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        records = [SeqRecord(Seq('A' * n), id=f'id_{n}') for n in (9, 12, 15, 18, 21)]
        groups = [runner.assign_records(records, 3, i) for i in range(3)]
        identifiers = [record.id for group in groups for record in group]
        self.assertCountEqual(identifiers, [record.id for record in records])
        self.assertEqual([r.id for r in groups[1]],
                         [r.id for r in runner.assign_records(list(reversed(records)), 3, 1)])

    def test_exclusive_lock_rejects_duplicates_and_can_be_reacquired(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'shard.lock'
            with runner.exclusive_lock(path):
                with self.assertRaises(OSError):
                    with runner.exclusive_lock(path, blocking=False):
                        pass
            with runner.exclusive_lock(path, blocking=False):
                pass


if __name__ == '__main__':
    unittest.main()
