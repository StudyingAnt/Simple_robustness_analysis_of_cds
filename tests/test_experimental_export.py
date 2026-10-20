"""Check scientific/audit binding, all-profile export and atomic campaign output."""
import argparse
from contextlib import ExitStack, redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

import numpy as np
from Bio import SeqIO

import test_bat_preparation as fixtures

EXPERIMENT = Path(__file__).resolve().parents[1] / '06_Bat_genome/experimental_apobec'
sys.path.insert(0, str(EXPERIMENT))
import export_results as exporter
import prepare_profiles as profiles
import run_experimental_apobec as runner
import aggregate_bat_results as bat_exporter
from validate_results import validate_run


class ExperimentalExportTests(unittest.TestCase):
    def setUp(self):
        fixture = fixtures.BatPreparationTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        fixture.fixture(genome='>chrP\nCATGGGGGTAGTAAG\n>chrN\nGTTACCCCTTTCATC\n')
        self.root = fixture.root
        self.cohorts = fixture.manifest
        body = json.loads(fixture.manifest.read_text())
        body['cohorts']['myoMyo6'] = body['cohorts']['molMol2'].copy()
        fixture.manifest.write_text(json.dumps(body))
        for original in (fixture.bed, fixture.genome):
            shutil.copy2(original, original.with_name(original.name.replace('molMol2', 'myoMyo6')))
        fixture.sources.write_text(json.dumps(dict(assets=[dict(filename=p.name, expected_bytes=p.stat().st_size,
            sha256=profiles.digest(p)) for p in fixture.raw.iterdir()])))
        self.campaign = dict(profiles_dir=str(self.root / 'profiles'), datasets={},
                             profile_ids=[row[1] for row in profiles.ENZYMES])
        identity = dict(profiles=[dict(id=r[1], sample_code=r[0], enzyme_species=r[2]) for r in profiles.ENZYMES])
        folder = self.root / 'profiles'
        folder.mkdir()
        (folder / 'profiles.json').write_text(json.dumps(identity))
        (folder / 'preparation_manifest.json').write_text('{}')
        self.campaign.update(profiles_sha256=profiles.digest(folder / 'profiles.json'),
            profile_manifest_sha256=profiles.digest(folder / 'preparation_manifest.json'))
        self.identity = identity
        raw_read_inputs = exporter.read_inputs
        self.reader = lambda processed, species, input_path=None, allow_subset=False, cohort_manifest=None: raw_read_inputs(
            processed, species, input_path, allow_subset, self.cohorts)
        for dataset in exporter.DATASETS:
            if dataset == 'human':
                path = self.root / 'human_inputs/human.fa'
                path.parent.mkdir()
                path.write_text('>ENST1|ENSG1|G|T|GENE-201|GENE|\nAATGAAATAAT\n'
                                '>ENST2|ENSG1|G|T|GENE-202|GENE|\nAATGTAATAAT\n')
                row = dict(input_path=str(path), input_sha256=profiles.digest(path), transcripts=2)
            else:
                processed = self.root / 'processed' / dataset
                with redirect_stdout(io.StringIO()):
                    fixtures.extractor.extract(fixture.raw, processed, dataset, fixture.sources, self.cohorts)
                inputs = self.reader(processed, dataset)
                row = dict(input_path=str(inputs['input_path']), input_sha256=inputs['input_sha256'], transcripts=2,
                           extraction_manifest_path=str(processed / 'extraction_manifest.json'))
            self.campaign['datasets'][dataset] = row
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        (self.inputs / 'inputs_manifest.json').write_text(json.dumps(self.campaign))
        self.runs, self.audits, self.output = [self.root / name for name in ('runs', 'audits', 'tables')]
        audit = dict(stage='audit', stage_passed=True, n_trajectories=3, datasets={},
                     campaign_sha256=profiles.digest(self.inputs / 'inputs_manifest.json'))
        for d, dataset in enumerate(exporter.DATASETS):
            run = self.runs / dataset
            (run / 'transcripts').mkdir(parents=True)
            (run / 'shards').mkdir()
            row = self.campaign['datasets'][dataset]
            config = dict(runner='rums-v3', mutation_model='experimental_apobec_NNC_v1', dataset=dataset,
                input_sha256=row['input_sha256'], signatures=self.campaign['profile_ids'], n_trajectories=3,
                profiles_sha256=self.campaign['profiles_sha256'],
                profile_preparation_sha256=self.campaign['profile_manifest_sha256'],
                criterion='protein_identity', strand='both', strand_bias=.5, cds_start=1,
                cds_end='len(sequence)-1', translation_table=1, internal_stops_retained_as_preprocessed=True,
                cosmic_version=None, genome_ref=None, shards=1, max_steps=100, master_seed=42,
                runner_python_sha256='same-source', package_source_sha256='same-package',
                profile_loader_sha256='same-loader', python='3.9.15', numpy='1.21.6', biopython='1.79',
                sinabro='3.0.0', assignment='greedy_length_balance_v1')
            fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
            (run / 'run_manifest.json').write_text(json.dumps(dict(configuration=config,
                config_sha256=fingerprint, input_path=row['input_path'])))
            for i, record in enumerate(SeqIO.parse(row['input_path'], 'fasta')):
                values = dict(counts=np.arange(18, dtype=np.int32).reshape(3, 6) + i + d,
                              codes=np.zeros((3, 6), dtype=np.uint8), seconds=np.ones(6))
                runner.save_checkpoint(run / 'transcripts' / (record.id.split('|')[0] + '.npz'),
                                       record.id, values, config, fingerprint)
            (run / 'shards/000.done.json').write_text(json.dumps(dict(config_sha256=fingerprint,
                assigned=2, skipped=0, newly_processed=2, incomplete_paths_in_new_results=0)))
            with redirect_stdout(io.StringIO()):
                self.assertTrue(validate_run(run, row['input_path'], self.audits / dataset)['validation_passed'])
            audit['datasets'][dataset] = dict(passed=True, transcripts=2, paths=36,
                report_sha256=profiles.digest(self.audits / dataset / 'validation_summary.json'))
        self.audit_path = self.audits / 'experimental_audit_manifest.json'
        self.audit_path.write_text(json.dumps(audit))
        self.args = argparse.Namespace(inputs_dir=self.inputs, run_root=self.runs, audit_root=self.audits,
                                       out_root=self.output, workers=1, n_traj=3, allow_subset=False)

    def export(self):
        # Only substitute the tiny frozen reference/source identities for test fixtures.
        # Real archive writers, audit readers and CSV exporters remain in use.
        with ExitStack() as stack:
            stack.enter_context(patch.object(exporter, 'load_inputs', return_value=self.campaign))
            stack.enter_context(patch.object(exporter, 'load_profiles', return_value=self.identity))
            stack.enter_context(patch.object(exporter, 'read_inputs', self.reader))
            stack.enter_context(patch.object(bat_exporter, 'read_inputs', self.reader))
            stack.enter_context(redirect_stdout(io.StringIO()))
            return exporter.export(self.args)

    def test_all_eighteen_summaries_preserve_every_profile_and_source_archive(self):
        before = {str(p): profiles.digest(p) for parent in (self.runs, self.audits)
                  for p in parent.rglob('*') if p.is_file()}
        result = self.export()
        self.assertTrue(result['stage_passed'])
        self.assertEqual((result['total_transcripts'], result['total_trajectories'], result['summary_rows']), (6, 108, 18))
        with (self.output / 'profile_comparison_summary.csv').open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 18)
        self.assertEqual([r['Signature'] for r in rows[:6]], self.campaign['profile_ids'])
        self.assertEqual(rows[3]['Signature'], 'HLmyoMyo6_APOBEC3A3')
        self.assertEqual([rows[i]['relation'] for i in (0, 6, 9, 12, 15)],
                         ['human_control', 'same_species', 'other_species', 'other_species', 'same_species'])
        self.assertEqual([float(r['mean_rums_across_transcripts']) for r in rows[:6]],
                         [6.5, 7.5, 8.5, 9.5, 10.5, 11.5])
        self.assertEqual(result['datasets']['human']['internal_stop_transcripts_retained'], 1)
        for dataset in ('molMol2', 'myoMyo6'):
            manifest = json.loads((self.output / dataset / 'aggregation_manifest.json').read_text())
            self.assertFalse(manifest['gene_mapping_available'])
            self.assertEqual(manifest['internal_stop_transcripts_retained'], 1)
        self.assertEqual(before, {str(p): profiles.digest(p) for parent in (self.runs, self.audits)
                                for p in parent.rglob('*') if p.is_file()})
        for entry in result['outputs']:
            self.assertEqual(entry['sha256'], profiles.digest(self.output / entry['filename']))

    def test_corruption_in_last_dataset_cannot_publish_partial_campaign(self):
        path = next((self.runs / 'myoMyo6/transcripts').glob('*.npz'))
        path.write_bytes(path.read_bytes() + b'changed since audit')
        with self.assertRaisesRegex(ValueError, 'changed since validation'):
            self.export()
        self.assertFalse(self.output.exists())

    def test_failed_or_unrelated_campaign_audit_is_rejected(self):
        original = json.loads(self.audit_path.read_text())
        for changes in (dict(stage_passed=False), dict(campaign_sha256='different'), dict(n_trajectories=1000)):
            self.audit_path.write_text(json.dumps(dict(original, **changes)))
            with self.assertRaisesRegex(ValueError, 'passing audit'):
                self.export()
            self.assertFalse(self.output.exists())

    def test_changed_audit_summary_is_rejected(self):
        path = self.audits / 'molMol2/validation_summary.json'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'audit summary hash'):
            self.export()
        self.assertFalse(self.output.exists())

    def test_wrong_profiles_or_scientific_criterion_are_rejected(self):
        path = self.runs / 'human/run_manifest.json'
        original = json.loads(path.read_text())
        for key, value in [('profiles_sha256', 'other'), ('criterion', 'blosum'),
                           ('signatures', self.campaign['profile_ids'][::-1])]:
            body = dict(original, configuration=dict(original['configuration'], **{key: value}))
            path.write_text(json.dumps(body))
            with self.assertRaises(ValueError):
                self.export()
            self.assertFalse(self.output.exists())

    def test_different_species_seed_is_rejected(self):
        path = self.runs / 'myoMyo6/run_manifest.json'
        body = json.loads(path.read_text())
        body['configuration']['master_seed'] = 43
        path.write_text(json.dumps(body))
        with self.assertRaisesRegex(ValueError, 'different simulation'):
            self.export()
        self.assertFalse(self.output.exists())

    def test_existing_output_and_nested_outputs_are_rejected(self):
        self.export()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.export()
        self.args.out_root = self.runs / 'nested'
        with self.assertRaisesRegex(ValueError, 'outside source'):
            self.export()


if __name__ == '__main__':
    unittest.main()
