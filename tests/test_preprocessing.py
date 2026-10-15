import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SCRIPT = Path(__file__).resolve().parents[1] / '01_data_downloading_and_preprocessing/extract_cds_plus.py'
spec = importlib.util.spec_from_file_location('preprocessing', SCRIPT)
preprocessing = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preprocessing)


def record(name, sequence, start=2, end=10):
    return SeqRecord(Seq(sequence), id=f'{name}|GENE|CDS:{start}-{end}', description='')


class PreprocessingTests(unittest.TestCase):
    def test_historical_filter_masks(self):
        cases = [
            (record('valid', 'CATGGCCTAAG'), 0),
            (record('length', 'CATGACTAAG', end=9), 1),
            (record('start', 'CCTGGCCTAAG'), 2),
            (record('stop', 'CATGGCCTACG'), 4),
            (record('id_PAR_Y', 'CATGGCCTAAG'), 8),
            (record('combined', 'CCTGGCTAAG', end=9), 3),
        ]
        for item, expected in cases:
            with self.subTest(id=item.id):
                self.assertEqual(preprocessing.classify(item)[0], expected)

    def test_missing_or_ambiguous_annotation_fails(self):
        for identifier in ['id|GENE', 'id|CDS:2-10|CDS:2-10', 'id|CDS:abc']:
            with self.subTest(id=identifier), self.assertRaises(ValueError):
                preprocessing.classify(SeqRecord(Seq('CATGGCCTAAG'), id=identifier))

    def test_coordinates_outside_transcript_fail(self):
        for start, end in [(0, 10), (11, 2), (2, 12)]:
            with self.subTest(start=start, end=end), self.assertRaises(ValueError):
                preprocessing.classify(record('bad', 'CATGGCCTAAG', start, end))

    def test_flanks_and_rejected_record_accounting(self):
        items = [
            record('valid', 'CATGGCCTAAG'),
            record('start', 'CCTGGCCTAAG'),
            record('id_PAR_Y', 'CATGGCCTAAG'),
            record('no_left_flank', 'ATGGCCTAAG', start=1, end=9),
            record('no_right_flank', 'CATGGCCTAA'),
        ]
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.fa'
            output = Path(directory) / 'processed'
            SeqIO.write(items, source, 'fasta')
            report = preprocessing.preprocess(source, output, 40)
            accepted = list(SeqIO.parse(output / 'gencode.v40.pc_transcripts.nopary.cdsplus.fa', 'fasta'))
            self.assertEqual(report['input_records'], 5)
            self.assertEqual(report['accepted_records'], 1)
            self.assertEqual(report['rejected_by_original_filters'], 2)
            self.assertEqual(report['otherwise_valid_but_missing_flank'], 2)
            self.assertEqual(str(accepted[0].seq), 'CATGGCCTAAG')
            self.assertEqual(str(accepted[0].seq[1:-1].translate()), 'MA*')
            self.assertEqual(sum(entry['fasta_records'] for entry in report['outputs']), 5)
            self.assertEqual(json.loads((output / 'preprocessing_manifest.json').read_text()), report)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.fa'
            output = Path(directory) / 'processed'
            output.mkdir()
            keep = output / 'keep.txt'
            keep.write_text('keep')
            SeqIO.write([record('valid', 'CATGGCCTAAG')], source, 'fasta')
            with self.assertRaises(FileExistsError):
                preprocessing.preprocess(source, output, 40)
            self.assertEqual(keep.read_text(), 'keep')

    def test_empty_input_does_not_publish_partial_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.fa'
            source.write_text('')
            output = Path(directory) / 'processed'
            with self.assertRaises(ValueError):
                preprocessing.preprocess(source, output, 40)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
