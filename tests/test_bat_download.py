"""Raw download checks without network requests or genome simulations."""
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / '06_Bat_genome/download_bat_inputs.py'
spec = importlib.util.spec_from_file_location('bat_download_target', SCRIPT)
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


class BatDownloadTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.config = self.root / 'sources.json'
        self.out = self.root / 'output'

    def source(self, raw, name='genome.fa.gz', expected_sha=None):
        asset = dict(filename=name, expected_bytes=len(raw), role='genome', url='https://example.invalid/data')
        if expected_sha:
            asset['sha256'] = expected_sha
        self.config.write_text(json.dumps(dict(assets=[asset])))

    def run_download(self, raw):
        with patch.object(downloader.urllib.request, 'urlopen', return_value=io.BytesIO(raw)):
            return downloader.download(self.out, self.config)

    def test_complete_gzip_download_checks_crc_and_records_hash(self):
        raw = gzip.compress(b'>contig\nATGTAAT\n')
        self.source(raw, expected_sha=hashlib.sha256(raw).hexdigest())
        report = self.run_download(raw)
        self.assertTrue(report['download_passed'])
        self.assertTrue(report['outputs'][0]['gzip_integrity_checked'])
        self.assertTrue(report['outputs'][0]['reference_sha256_checked'])
        self.assertEqual((self.out / 'genome.fa.gz').read_bytes(), raw)
        with self.assertRaises(FileExistsError):
            self.run_download(raw)

    def test_hash_mismatch_cannot_publish_output(self):
        raw = b'annotation'
        self.source(raw, name='genes.bed', expected_sha='0'*64)
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.run_download(raw)
        self.assertFalse(self.out.exists())

    def test_size_mismatch_cannot_publish_output(self):
        self.source(b'expected bytes', name='genes.bed')
        with self.assertRaisesRegex(ValueError, 'Size mismatch'):
            self.run_download(b'short')
        self.assertFalse(self.out.exists())

    def test_corrupt_gzip_cannot_publish_output(self):
        raw = b'not a gzip file'
        self.source(raw)
        with self.assertRaises(OSError):
            self.run_download(raw)
        self.assertFalse(self.out.exists())


if __name__ == '__main__':
    unittest.main()
