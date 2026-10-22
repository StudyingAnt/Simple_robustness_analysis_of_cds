"""Experimental profile provenance, normalization, checkpoint and allocation checks."""
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import zipfile

import numpy as np

EXPERIMENT = Path(__file__).resolve().parents[1] / '06_Bat_genome/experimental_apobec'
sys.path.insert(0, str(EXPERIMENT))
import prepare_profiles as profiles
import run_campaign as campaign
import run_experimental_apobec as runner
sys.path.insert(0, str(EXPERIMENT / 'author_hpc'))
import validate_experimental_pilot as pilot
from validate_results import audit_checkpoint

XMLNS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'


def worksheet(values, formula=False):
    root = ET.Element('worksheet', xmlns=XMLNS)
    data = ET.SubElement(root, 'sheetData')
    for key, value in values.items():
        cell = ET.SubElement(data, 'c', r=key)
        if isinstance(value, str):
            cell.set('t', 'inlineStr')
            ET.SubElement(ET.SubElement(cell, 'is'), 't').text = value
        else:
            ET.SubElement(cell, 'v').text = str(value)
    if formula:
        ET.SubElement(ET.SubElement(data, 'c', r='Z1'), 'f').text = '1+1'
    return ET.tostring(root)


class ExperimentalProfileTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.workbook = self.root / 'input.xlsx'
        self.out = self.root / 'profiles'

    def fixture(self, change=None, formula=False):
        data = {'C23': 'C>T counts and Percentage are average for each motif'}
        key = {}
        for i, (sample, gene, _, count_col, total_col, percent_col) in enumerate(profiles.ENZYMES):
            key[f'A{i+3}'], key[f'B{i+3}'] = sample, gene
            data[f'{total_col}2'] = sample
            for j, motif in enumerate(profiles.MOTIFS):
                row = j + 4
                data[f'B{row}'], data[f'C{row}'] = motif, j % 3 + 1
                count, total = (i + 1) * (j + 1), 10000
                data[f'{count_col}{row}'], data[f'{total_col}{row}'] = count, total
                data[f'{percent_col}{row}'] = 100 * count / total
        if change:
            change(data, key)
        with zipfile.ZipFile(self.workbook, 'w') as archive:
            archive.writestr('xl/worksheets/sheet1.xml', worksheet(data, formula))
            archive.writestr('xl/worksheets/sheet2.xml', worksheet(key))
        return profiles.digest(self.workbook)

    def prepare(self, **kwargs):
        return profiles.prepare(self.workbook, self.out, expected_sha=profiles.digest(self.workbook), **kwargs)

    def test_means_become_per_site_weights_without_substrate_count_multiplication(self):
        sha = self.fixture()
        self.prepare()
        with patch.object(profiles, 'WORKBOOK_SHA256', sha):
            body = profiles.load_profiles(self.out)
        self.assertEqual(len(body['profiles']), 6)
        self.assertEqual(body['profiles'][3]['id'], 'HLmyoMyo6_APOBEC3A3')
        np.testing.assert_allclose(body['profiles'][0]['weights'], np.arange(1, 17) / 136.)
        self.assertFalse(body['substrate_motif_counts_multiplied'])
        self.assertFalse(body['historical_probability_csv_recovered'])

    def test_wrong_original_checksum_is_rejected_without_publishing(self):
        self.fixture()
        with self.assertRaisesRegex(ValueError, 'workbook SHA256'):
            profiles.prepare(self.workbook, self.out)
        self.assertFalse(self.out.exists())

    def test_formulas_are_not_executed_or_accepted(self):
        self.fixture(formula=True)
        with self.assertRaisesRegex(ValueError, 'formulas'):
            self.prepare()

    def test_sample_identity_mismatch_is_rejected(self):
        self.fixture(change=lambda _, key: key.update(B3='other_enzyme'))
        with self.assertRaisesRegex(ValueError, 'identities'):
            self.prepare()

    def test_inconsistent_percentage_is_rejected(self):
        self.fixture(change=lambda data, _: data.update(F4=42))
        with self.assertRaisesRegex(ValueError, 'Percentage differs'):
            self.prepare()

    def test_duplicate_motif_is_rejected(self):
        self.fixture(change=lambda data, _: data.update(B5=data['B4']))
        with self.assertRaisesRegex(ValueError, 'Missing/duplicate'):
            self.prepare()

    def test_outputs_are_immutable_and_changed_weights_cannot_pass_hash_check(self):
        sha = self.fixture()
        self.prepare()
        with self.assertRaisesRegex(ValueError, 'new output directory'):
            self.prepare()
        path = self.out / 'profiles.json'
        path.write_text(path.read_text() + '\n')
        with patch.object(profiles, 'WORKBOOK_SHA256', sha), self.assertRaisesRegex(ValueError, 'profile hash'):
            profiles.load_profiles(self.out)


class ExperimentalCheckpointTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.config = dict(runner='rums-v3', signatures=[r[1] for r in profiles.ENZYMES],
                           n_trajectories=3, max_steps=100)
        self.fingerprint = hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()
        self.path = self.root / 'transcript.npz'
        self.values = dict(counts=np.arange(18, dtype=np.int32).reshape(3, 6),
                           codes=np.zeros((3, 6), dtype=np.uint8), seconds=np.ones(6))

    def test_saved_statistics_pass_independent_archive_audit_and_refuse_overwrite(self):
        self.assertEqual(runner.save_checkpoint(self.path, 'transcript', self.values, self.config, self.fingerprint), 0)
        result = audit_checkpoint((self.path, 'transcript'), self.config, self.fingerprint)
        self.assertTrue(result['structurally_valid'])
        self.assertEqual(result['complete_paths'], 18)
        self.assertTrue(runner.check_existing(self.path, 'transcript', self.config, self.fingerprint))
        with self.assertRaises(FileExistsError):
            runner.save_checkpoint(self.path, 'transcript', self.values, self.config, self.fingerprint)

    def test_incomplete_path_has_nan_statistics_and_is_reported_in_audit(self):
        self.values['counts'][0, 0] = -1
        self.values['codes'][0, 0] = 1
        self.assertEqual(runner.save_checkpoint(self.path, 'transcript', self.values, self.config, self.fingerprint), 1)
        result = audit_checkpoint((self.path, 'transcript'), self.config, self.fingerprint)
        self.assertTrue(result['structurally_valid'])
        self.assertEqual(result['no_available_event_paths'], 1)
        with np.load(self.path, allow_pickle=False) as saved:
            self.assertTrue(np.isnan(saved['mean'][0]))

    def test_changed_existing_checkpoint_is_rejected(self):
        runner.save_checkpoint(self.path, 'transcript', self.values, self.config, self.fingerprint)
        with self.assertRaisesRegex(ValueError, 'Invalid existing checkpoint'):
            runner.check_existing(self.path, 'another', self.config, self.fingerprint)

    def test_pilot_comparison_checks_each_path_not_only_means(self):
        left, right = self.root / 'left', self.root / 'right'
        for directory in (left, right):
            (directory / 'transcripts').mkdir(parents=True)
            runner.save_checkpoint(directory / 'transcripts/transcript.npz', 'transcript', self.values,
                                   self.config, self.fingerprint)
        self.assertEqual(pilot.compare(left, right, 3, self.config['signatures'], expected_records=1), 18)
        (right / 'transcripts/transcript.npz').unlink()
        changed = dict(self.values, counts=self.values['counts'][::-1].copy())
        runner.save_checkpoint(right / 'transcripts/transcript.npz', 'transcript', changed, self.config, self.fingerprint)
        with self.assertRaisesRegex(ValueError, 'Per-path'):
            pilot.compare(left, right, 3, self.config['signatures'], expected_records=1)


class ExperimentalAllocationTests(unittest.TestCase):
    def test_one_fifty_task_submission_covers_all_three_hundred_cases_once(self):
        plan = [campaign.task_cases(i, 50, 100) for i in range(50)]
        self.assertTrue(all(len(row) == 6 for row in plan))
        cases = [case for row in plan for case in row]
        self.assertEqual(len(cases), len(set(cases)))
        self.assertEqual(set(cases), {(dataset, shard) for dataset in campaign.DATASETS for shard in range(100)})
        self.assertEqual(plan[0], [(dataset, shard) for shard in (0, 50) for dataset in campaign.DATASETS])

    def test_invalid_allocation_parameters_fail(self):
        for args in [(-1, 50, 100), (50, 50, 100), (0, 0, 100), (0, 101, 100)]:
            with self.assertRaises(ValueError):
                campaign.task_cases(*args)


class ExperimentalPilotGateTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        (self.inputs / 'inputs_manifest.json').write_text('{}')
        self.pilot = self.root / 'pilot'
        (self.pilot / 'inputs').mkdir(parents=True)
        self.campaign = dict(profiles_sha256='profiles', profile_ids=['six-enzymes'])
        report = dict(validation_passed=True, n_trajectories=100, datasets={},
                      campaign_sha256=profiles.digest(self.inputs / 'inputs_manifest.json'), **self.campaign)
        for dataset in campaign.DATASETS:
            fasta = self.pilot / 'inputs' / f'{dataset}.fa'
            fasta.write_text('>test\nACG\n')
            report['datasets'][dataset] = dict(results_equal=True, all_paths_complete=True,
                observed_parallel_execution=True, audits=dict(serial=True, parallel=True),
                input_sha256=profiles.digest(fasta))
            for mode in ('serial', 'parallel'):
                path = self.pilot / dataset / mode
                path.mkdir(parents=True)
                config = dict(runner_python_sha256=profiles.digest(EXPERIMENT / 'run_experimental_apobec.py'),
                    profile_loader_sha256=profiles.digest(EXPERIMENT / 'prepare_profiles.py'),
                    package_source_sha256='package', master_seed=42, max_steps=10000)
                (path / 'run_manifest.json').write_text(json.dumps(dict(configuration=config)))
        (self.pilot / 'experimental_pilot_report.json').write_text(json.dumps(report))

    def test_current_code_package_and_settings_are_accepted(self):
        campaign.check_pilot(self.pilot, self.inputs, self.campaign, 'package', 42, 10000)

    def test_changed_package_seed_or_cap_requires_new_pilot(self):
        for settings in [('other-package', 42, 10000), ('package', 43, 10000), ('package', 42, 100)]:
            with self.subTest(settings=settings), self.assertRaisesRegex(ValueError, 'fresh pilot'):
                campaign.check_pilot(self.pilot, self.inputs, self.campaign, *settings)

    def test_changed_runner_is_not_covered_by_previous_pilot(self):
        path = self.pilot / 'human' / 'parallel' / 'run_manifest.json'
        body = json.loads(path.read_text())
        body['configuration']['runner_python_sha256'] = 'older-runner'
        path.write_text(json.dumps(body))
        with self.assertRaisesRegex(ValueError, 'fresh pilot'):
            campaign.check_pilot(self.pilot, self.inputs, self.campaign, 'package', 42, 10000)


if __name__ == '__main__':
    unittest.main()
