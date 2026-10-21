"""Integrity, matched identity, and output provenance of recovered S11 estimates."""
import argparse
import csv
import contextlib
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

PAPER = Path(__file__).resolve().parents[1]/'07_manuscript_validation'
sys.path.insert(0, str(PAPER))
import analyze_s11_historical as s11


class HistoricalS11Tests(unittest.TestCase):
    def fixture(self, root, constant=False):
        raw = root/'raw'; raw.mkdir()
        fasta = raw/'human.fa'
        # Input order differs from historical row order on purpose.
        fasta.write_text('>tx2|g|x|x|B-201\nCATGGCTTAAG\n>tx1|g|x|x|A-201\nCATGGCCTAAG\n')
        stream = io.StringIO(newline='')
        writer = csv.writer(stream, lineterminator='\n'); writer.writerow(s11.FIELDS)
        writer.writerow([0, 'tx1', 'A-201', 0, 2])
        writer.writerow([1, 'tx2', 'B-201', 0 if constant else 1, 2 if constant else 3])
        data = stream.getvalue().encode()
        values = raw/'values.csv.gz'; values.write_bytes(gzip.compress(data, mtime=0))
        body = dict(schema_version=1, historical_results_only=True, signature='SBS2',
            transcripts=2, values_filename=values.name, values_sha256=s11.digest(values),
            uncompressed_csv_sha256=hashlib.sha256(data).hexdigest(),
            human_input_sha256=s11.digest(fasta), historical_script_n_sim=1000,
            historical_criterion=dict(matrix='BLOSUM80', threshold=20),
            historical_seeds_recovered=False, historical_plotted_populations=dict(identity=7, blosum=2))
        manifest = raw/'manifest.json'; manifest.write_text(json.dumps(body))
        return argparse.Namespace(input=fasta, manifest=manifest, values=None, out_dir=root/'result')

    def test_bundled_recovered_population_and_checksum(self):
        manifest = json.loads((PAPER/'source_data/s11_historical_manifest.json').read_text())
        rows = s11.read_values(PAPER/'source_data'/manifest['values_filename'], manifest)
        self.assertEqual(len(rows), 3740)
        self.assertEqual(rows[0]['transcript_name'], 'OR4F5-201')
        self.assertEqual(rows[-1]['transcript_name'], 'PLEKHA3-201')
        self.assertAlmostEqual(sum(r['historical_blosum_rums'] for r in rows)/len(rows), 3.8720010695187166)

    def test_matched_output_preserves_order_sources_and_historical_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp))
            before = {p:s11.digest(p) for p in args.manifest.parent.iterdir()}
            with contextlib.redirect_stdout(io.StringIO()):
                result = s11.run(args)
            with (args.out_dir/'s11_matched_transcript_values.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([r['transcript_id'] for r in rows], ['tx1', 'tx2'])
            self.assertEqual(result['identity']['mean'], .5)
            self.assertEqual(result['blosum']['mean'], 2.5)
            self.assertTrue(result['analysis_passed'])
            self.assertFalse(result['fresh_simulation'])
            self.assertFalse(result['corrected_package_blosum_recalculation'])
            self.assertEqual(result['input_sequences_verified'], 2)
            self.assertTrue((args.out_dir/'figureS11_historical_matched.pdf').exists())
            self.assertTrue(all(s11.digest(p)==value for p,value in before.items()))
            for entry in result['outputs']:
                self.assertEqual(s11.digest(args.out_dir/entry['filename']),entry['sha256'])

    def test_changed_values_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp))
            (args.manifest.parent/'values.csv.gz').write_bytes(gzip.compress(b'changed'))
            with self.assertRaisesRegex(ValueError, 'values checksum'):
                s11.run(args)
            self.assertFalse(args.out_dir.exists())

    def test_changed_reference_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp)); args.input.write_text('>wrong\nAAAA\n')
            with self.assertRaisesRegex(ValueError, 'FASTA checksum'):
                s11.run(args)
            self.assertFalse(args.out_dir.exists())

    def test_wrong_name_is_rejected_even_with_valid_fasta_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp)); args.input.write_text(args.input.read_text().replace('A-201','X-201'))
            body = json.loads(args.manifest.read_text()); body['human_input_sha256']=s11.digest(args.input)
            args.manifest.write_text(json.dumps(body))
            with self.assertRaisesRegex(ValueError, 'ID/name mismatch'):
                s11.run(args)

    def test_existing_output_and_nested_input_destinations_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp)); args.out_dir.mkdir()
            with self.assertRaisesRegex(ValueError, 'new S11 output'):
                s11.run(args)
            args.out_dir=args.input.parent/'nested'
            with self.assertRaisesRegex(ValueError, 'nested inside'):
                s11.run(args)

    def test_zero_and_constant_estimates_render_without_kde_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.fixture(Path(tmp), constant=True)
            with contextlib.redirect_stdout(io.StringIO()):
                report = s11.run(args)
            self.assertEqual(report['identity']['mean'], 0)
            self.assertEqual(report['blosum']['mean'], 2)


if __name__ == '__main__':
    unittest.main()
