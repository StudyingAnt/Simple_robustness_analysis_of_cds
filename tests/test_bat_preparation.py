"""Coordinate, cohort, failure-publication and pilot-integrity checks."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / '06_Bat_genome' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


extractor = load('bat_extract_test_target', 'extract_bat_cds.py')
pilot = load('bat_pilot_test_target', 'validate_bat_pilot.py')


class BatPreparationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.raw = self.root / 'raw'
        self.raw.mkdir()
        self.out = self.root / 'out'
        self.bed = self.raw / 'HLmolMol2.Bat1Kannotation.bed'
        self.genome = self.raw / 'HLmolMol2.fa.gz'
        self.sources = self.root / 'sources.json'
        self.cohort = self.root / 'cohort.txt'
        self.manifest = self.root / 'cohorts.json'

    def fixture(self, bed=None, genome=None, ids=('positive', 'negative')):
        self.bed.write_text(bed or
            'chrP\t1\t14\tpositive\t0\t+\t1\t14\t0\t2\t3,6,\t0,7,\n'
            'chrN\t1\t14\tnegative\t0\t-\t1\t14\t0\t2\t3,6\t0,7\n'
            'chrP\t1\t14\tunselected\t0\t+\t1\t14\t0\t2\t3,6,\t0,7,\n')
        with gzip.open(self.genome, 'wt') as stream:
            stream.write(genome or '>chrP\nCATGGGGGAAATAAG\n>chrN\nGTTACCCCTTTCATC\n')
        self.sources.write_text(json.dumps(dict(assets=[dict(filename=p.name, expected_bytes=p.stat().st_size,
            sha256=extractor.digest(p)) for p in (self.bed, self.genome)])))
        self.cohort.write_text('\n'.join(ids) + '\n')
        self.manifest.write_text(json.dumps(dict(cohorts=dict(molMol2=dict(
            filename=self.cohort.name, transcripts=len(ids), sha256=extractor.digest(self.cohort))))))

    def extract(self):
        return extractor.extract(self.raw, self.out, 'molMol2', self.sources, self.manifest)

    def test_spliced_both_strands_have_expected_oriented_flanks_and_frozen_ids(self):
        self.fixture()
        report = self.extract()
        self.assertTrue(report['extraction_passed'])
        actual = {r.id: str(r.seq) for r in SeqIO.parse(self.out / report['fasta_filename'], 'fasta')}
        self.assertEqual(actual, dict(positive='CATGAAATAAG', negative='GATGAAATAAC'))
        self.assertEqual(report['valid_transcripts'], 2)
        with self.assertRaises(FileExistsError):
            self.extract()

    def test_missing_cohort_id_is_rejected(self):
        self.fixture(ids=('absent',))
        with self.assertRaisesRegex(ValueError, 'absent from BED'):
            self.extract()
        self.assertFalse(self.out.exists())

    def test_changed_source_cannot_silently_change_inputs(self):
        self.fixture()
        self.bed.write_text(self.bed.read_text().replace('chrP', 'chrZ'))
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.extract()

    def test_changed_cohort_cannot_silently_change_population(self):
        self.fixture()
        self.cohort.write_text('unselected\n')
        with self.assertRaisesRegex(ValueError, 'Cohort SHA256'):
            self.extract()

    def test_cds_outer_bounds_mismatch_is_rejected(self):
        self.fixture(bed='chrP\t0\t14\tpositive\t0\t+\t1\t14\t0\t1\t14,\t0,\n', ids=('positive',))
        with self.assertRaisesRegex(ValueError, 'different transcript/CDS bounds'):
            self.extract()

    def test_unknown_base_publishes_failed_qc_without_final_fasta(self):
        self.fixture(genome='>chrP\nCATGGGGGAANTAAG\n>chrN\nGTTACCCCTTTCATC\n')
        with self.assertRaisesRegex(ValueError, 'invalid cohort transcripts'):
            self.extract()
        report = json.loads((self.out / 'extraction_manifest.json').read_text())
        self.assertFalse(report['extraction_passed'])
        self.assertEqual(report['errors'], dict(non_ACGT=1))
        self.assertTrue((self.out / 'transcript_qc.csv').exists())
        self.assertEqual(list(self.out.glob('*.fa')), [])

    def test_internal_stop_is_retained_and_flagged(self):
        self.fixture(genome='>chrP\nCATGGGGGTAGTAAG\n>chrN\nGTTACCCCTTTCATC\n')
        report = self.extract()
        self.assertTrue(report['extraction_passed'])
        self.assertEqual(report['internal_stop_transcripts'], 1)

    def test_missing_contig_publishes_failed_qc(self):
        self.fixture(genome='>chrP\nCATGGGGGAAATAAG\n')
        with self.assertRaisesRegex(ValueError, 'invalid cohort transcripts'):
            self.extract()
        report = json.loads((self.out / 'extraction_manifest.json').read_text())
        self.assertEqual(report['errors'], dict(missing_contig=1))

    def test_missing_genomic_flank_is_rejected(self):
        row = dict(transcript_id='edge', cds_start=0, cds_end=9, blocks=[(0, 9)], strand='+')
        with self.assertRaisesRegex(ValueError, 'Missing flank'):
            extractor.extract_sequence(Seq('ATGAAATAAG'), row)

    def test_corrupt_gzip_is_not_accepted_even_if_source_hash_matches(self):
        self.fixture()
        self.genome.write_bytes(b'not gzip')
        self.sources.write_text(json.dumps(dict(assets=[dict(filename=p.name, expected_bytes=p.stat().st_size,
            sha256=extractor.digest(p)) for p in (self.bed, self.genome)])))
        with self.assertRaises(OSError):
            self.extract()
        self.assertFalse(self.out.exists())

    def test_pilot_checks_fasta_hash_and_excludes_internal_stops_only_from_pilot(self):
        self.fixture(genome='>chrP\nCATGGGGGTAGTAAG\n>chrN\nGTTACCCCTTTCATC\n')
        self.extract()
        _, selected, population = pilot.select_records(self.out, 1)
        self.assertEqual([r.id for r in selected], ['negative'])
        self.assertEqual(population, dict(total=2, eligible=1, internal_stop_excluded=1))
        fasta = self.out / 'molMol2.cdsplus.fa'
        fasta.write_text(fasta.read_text().replace('ATG', 'ATA'))
        with self.assertRaisesRegex(ValueError, 'FASTA SHA256'):
            pilot.select_records(self.out, 1)


class PilotConsistencyTests(unittest.TestCase):
    def test_quantile_selection_is_order_independent_and_covers_lengths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / 'molMol2.cdsplus.fa'
            records = [SeqRecord(Seq('C' + 'ATG' + 'AAA' * i + 'TAA' + 'G'),
                                 id=f't{i}', description='') for i in range(24)]
            SeqIO.write(reversed(records), fasta, 'fasta')
            manifest = dict(extraction_passed=True, expected_transcripts=24, valid_transcripts=24,
                            fasta_filename=fasta.name, fasta_sha256=pilot.digest(fasta))
            (root / 'extraction_manifest.json').write_text(json.dumps(manifest))
            _, chosen, _ = pilot.select_records(root, 4)
            self.assertEqual([r.id for r in chosen], ['t3', 't9', 't15', 't21'])

    def test_per_path_comparison_rejects_changed_result(self):
        with tempfile.TemporaryDirectory() as temporary:
            roots = [Path(temporary) / label for label in ('serial', 'parallel')]
            for root in roots:
                (root / 'transcripts').mkdir(parents=True)
                np.savez(root / 'transcripts/id.npz', accepted_steps=np.ones((2, 3), dtype='int32'),
                         terminal_codes=np.zeros((2, 3), dtype='uint8'),
                         metadata=json.dumps(dict(config_sha256='same', signatures=pilot.SIGNATURES)))
            self.assertEqual(pilot.compare(*roots, 1, 2), 6)
            np.savez(roots[1] / 'transcripts/id.npz', accepted_steps=np.zeros((2, 3), dtype='int32'),
                     terminal_codes=np.zeros((2, 3), dtype='uint8'),
                     metadata=json.dumps(dict(config_sha256='same', signatures=pilot.SIGNATURES)))
            with self.assertRaisesRegex(ValueError, 'Per-path'):
                pilot.compare(*roots, 1, 2)


if __name__ == '__main__':
    unittest.main()
