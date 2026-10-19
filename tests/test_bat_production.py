"""Production array mapping, provenance and bat-specific audited table export."""
import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
from Bio import SeqIO

import test_bat_preparation as preparation

BAT = Path(__file__).resolve().parents[1] / '06_Bat_genome'
sys.path.insert(0, str(BAT))
import bat_inputs
import aggregate_bat_results as exporter
import run_bat_cosmic as driver
import run_rums as runner
import validate_results as validator


class ArrayTests(unittest.TestCase):
    def test_array_covers_each_species_shard_once(self):
        mapped = [bat_inputs.array_case(i, 100) for i in range(200)]
        self.assertEqual(len(set(mapped)), 200)
        self.assertEqual(mapped[:4], [('molMol2', 0), ('myoMyo6', 0), ('molMol2', 1), ('myoMyo6', 1)])
        self.assertEqual(set(mapped), {(s, i) for s in bat_inputs.SPECIES for i in range(100)})

    def test_invalid_array_range_fails(self):
        for index, shards in [(-1, 100), (200, 100), (0, 0), (0, -1)]:
            with self.assertRaises(ValueError):
                bat_inputs.array_case(index, shards)


class BatProductionTests(unittest.TestCase):
    def setUp(self):
        # Reuse the coordinate fixture; no network or biological simulation is needed.
        fixture = preparation.BatPreparationTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        fixture.fixture(genome='>chrP\nCATGGGGGTAGTAAG\n>chrN\nGTTACCCCTTTCATC\n')
        with redirect_stdout(io.StringIO()):
            fixture.extract()
        self.fixture = fixture
        self.base = fixture.root
        self.run = self.base / 'run'
        self.audit = self.base / 'audit'
        self.tables = self.base / 'tables'
        self.make_run()

    def inputs(self, input_path=None, allow_subset=False):
        return bat_inputs.read_inputs(self.fixture.out, 'molMol2', input_path, allow_subset, self.fixture.manifest)

    def make_run(self):
        inputs = self.inputs()
        self.fasta = inputs['input_path']
        self.config = dict(runner='rums-v3', input_sha256=inputs['input_sha256'],
            n_trajectories=3, signatures=bat_inputs.SIGNATURES, shards=1, max_steps=100,
            assignment='greedy_length_balance_v1', criterion='protein_identity', cosmic_version='3.3.1',
            genome_ref='GRCh37', cds_start=1, cds_end='len(sequence)-1', translation_table=1,
            internal_stops_retained_as_preprocessed=True)
        self.fingerprint = hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()
        (self.run / 'transcripts').mkdir(parents=True)
        (self.run / 'shards').mkdir()
        (self.run / 'run_manifest.json').write_text(json.dumps(dict(config_sha256=self.fingerprint,
            configuration=self.config, input_path=str(self.fasta))))
        with patch.multiple(runner, OUTPUT=self.run, N_TRAJ=3, SIGNATURES=bat_inputs.SIGNATURES,
                            FINGERPRINT=self.fingerprint, create=True):
            for i, row in enumerate(inputs['annotations']):
                runner.save_transcript(row['transcript_id'], row['transcript_id'], dict(
                    counts=np.arange(24, dtype=np.int32).reshape(3, 8) + i,
                    codes=np.zeros((3, 8), dtype=np.uint8), seconds=np.ones(8)))
        (self.run / 'shards/000.done.json').write_text(json.dumps(dict(config_sha256=self.fingerprint,
            assigned=2, skipped=0, newly_processed=2, incomplete_paths_in_new_results=0)))
        with redirect_stdout(io.StringIO()):
            self.assertTrue(validator.validate_run(self.run, self.fasta, self.audit)['validation_passed'])

    def export(self, **kwargs):
        with redirect_stdout(io.StringIO()):
            return exporter.aggregate_bat(self.run, self.audit, self.fixture.out, self.tables, 'molMol2',
                n_traj=3, cohort_manifest=self.fixture.manifest, **kwargs)

    def test_frozen_input_has_species_coordinates_and_no_invented_gene_mapping(self):
        inputs = self.inputs()
        self.assertEqual(inputs['scope'], 'full_historical_cohort')
        self.assertEqual({r['species'] for r in inputs['annotations']}, {'molMol2'})
        self.assertNotIn('gene_id', inputs['annotations'][0])
        self.assertEqual(len(inputs['annotations']), 2)

    def test_qc_changes_are_rejected(self):
        path = self.fixture.out / 'transcript_qc.csv'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'QC hash'):
            self.inputs()

    def test_subset_requires_explicit_flag_and_exact_sequence_identity(self):
        subset = self.base / 'subset.fa'
        records = list(SeqIO.parse(self.fasta, 'fasta'))
        SeqIO.write(records[:1], subset, 'fasta')
        with self.assertRaisesRegex(ValueError, 'allow-subset'):
            self.inputs(subset)
        inputs = self.inputs(subset, True)
        self.assertEqual(inputs['scope'], 'explicit_subset_of_historical_cohort')
        self.assertEqual(len(inputs['annotations']), 1)
        subset.write_text(subset.read_text().replace('ATG', 'ATA'))
        with self.assertRaisesRegex(ValueError, 'Subset sequence differs'):
            self.inputs(subset, True)

    def test_wrong_scientific_configuration_is_rejected(self):
        for key, value in [('signatures', bat_inputs.SIGNATURES[:-1]), ('genome_ref', 'GRCh38'),
                           ('criterion', 'blosum'), ('cds_start', 0)]:
            with self.assertRaisesRegex(ValueError, 'scientific configuration'):
                bat_inputs.check_config(dict(self.config, **{key: value}), self.inputs(), 3)
        with self.assertRaisesRegex(ValueError, 'number of trajectories'):
            bat_inputs.check_config(self.config, self.inputs(), 1000)

    def test_export_keeps_bat_ids_all_eight_profiles_stops_and_bound_hashes(self):
        result = self.export()
        self.assertTrue(result['aggregation_passed'])
        self.assertEqual(result['population'], 'full_historical_cohort')
        self.assertEqual(result['total_trajectories'], 48)
        self.assertEqual(result['internal_stop_transcripts_retained'], 1)
        self.assertFalse(result['gene_mapping_available'])
        self.assertNotIn('unique_genes', result)
        self.assertEqual((self.tables / 'rums_mean_by_transcript_id.csv').read_text(),
                         (self.tables / 'molMol2_rums_profile.csv').read_text())
        header = (self.tables / 'molMol2_rums_profile.csv').read_text().splitlines()[0]
        self.assertEqual(header.split(','), ['Transcript', *bat_inputs.SIGNATURES])
        for entry in result['outputs']:
            self.assertEqual(entry['sha256'], bat_inputs.digest(self.tables / entry['filename']))

    def test_export_rejects_changed_archive_after_audit(self):
        path = self.run / 'transcripts/positive.npz'
        path.write_bytes(path.read_bytes() + b'changed')
        with self.assertRaisesRegex(ValueError, 'changed since validation'):
            self.export()
        self.assertFalse(self.tables.exists())

    def test_export_rejects_wrong_audit_manifest_or_inventory(self):
        path = self.audit / 'validation_summary.json'
        report = json.loads(path.read_text())
        report['manifest_sha256'] = '0' * 64
        path.write_text(json.dumps(report))
        with self.assertRaisesRegex(ValueError, 'Audit belongs'):
            self.export()

    def test_export_rejects_failed_audit(self):
        path = self.audit / 'validation_summary.json'
        report = json.loads(path.read_text())
        report['validation_passed'] = False
        path.write_text(json.dumps(report))
        with self.assertRaisesRegex(ValueError, 'passing complete archive audit'):
            self.export()

    def test_export_refuses_overwrite_and_nested_outputs(self):
        self.export()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.export()
        with self.assertRaisesRegex(ValueError, 'outside'):
            exporter.aggregate_bat(self.run, self.audit, self.fixture.out, self.run / 'nested', 'molMol2')

    def test_driver_uses_existing_runner_and_refuses_changed_resume_settings(self):
        pilot = self.base / 'pilot.json'
        pilot.write_text('{}')
        args = argparse.Namespace(workers=4, n_traj=1000, seed=42, max_steps=10000,
            array_index=2, shards_per_species=100, processed_root=self.base / 'processed',
            pilot_dir=self.base / 'pilot', out_root=self.base / 'production', sinabro_source=None)
        with patch.object(driver, 'read_inputs', return_value=self.inputs()), \
             patch.object(driver, 'check_pilot', return_value=pilot), \
             patch.object(driver.subprocess, 'run') as launched, redirect_stdout(io.StringIO()):
            driver.run(args)
            command = launched.call_args.args[0]
            self.assertEqual(command[command.index('--shard-index') + 1], '1')
            self.assertEqual(command[command.index('--expected-records') + 1], '2')
            self.assertEqual(command[command.index('--signatures') + 1:], bat_inputs.SIGNATURES)
            driver.run(args)
            args.n_traj = 100
            with self.assertRaisesRegex(ValueError, 'provenance/configuration changed'):
                driver.run(args)


class PilotGateTests(unittest.TestCase):
    def test_pilot_must_pass_and_belong_to_current_extraction(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            pilot, processed = base / 'pilot', base / 'processed'
            report = dict(validation_passed=True, n_trajectories=100, signatures=['SBS2', 'SBS4', 'SBS13'], species={})
            for species in bat_inputs.SPECIES:
                (processed / species).mkdir(parents=True)
                (pilot / species).mkdir(parents=True)
                (processed / species / 'extraction_manifest.json').write_text('{}')
                (pilot / species / 'selected_transcripts.fa').write_text('>id\nAATGTAAT\n')
                report['species'][species] = dict(results_equal=True, all_paths_complete=True,
                    observed_parallel_execution=True, audits=dict(serial=True, parallel=True),
                    extraction_manifest_sha256=bat_inputs.digest(processed / species / 'extraction_manifest.json'),
                    selected_fasta_sha256=bat_inputs.digest(pilot / species / 'selected_transcripts.fa'))
            path = pilot / 'pilot_report.json'
            path.write_text(json.dumps(report))
            self.assertEqual(bat_inputs.check_pilot(pilot, processed), path)
            report['n_trajectories'] = 20
            path.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'at least 100'):
                bat_inputs.check_pilot(pilot, processed)
            report['n_trajectories'] = 100
            path.write_text(json.dumps(report))
            (processed / 'myoMyo6/extraction_manifest.json').write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError, 'different extraction manifest'):
                bat_inputs.check_pilot(pilot, processed)


if __name__ == '__main__':
    unittest.main()
