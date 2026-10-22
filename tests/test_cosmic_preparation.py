"""Audit TAR/TSV handling without redistributing COSMIC observations."""
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "01_data_downloading_and_preprocessing"))
import prepare_cosmic_brca1 as prep


class CosmicPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.archive = self.root / "Cosmic_MutantCensus_Tsv_v101_GRCh38.tar"
        self.manifest = self.root / "manifest.json"
        self.out = self.root / "processed"
        self.columns = ["GENE_SYMBOL", "TRANSCRIPT_ACCESSION", "COSMIC_PHENOTYPE_ID", "MUTATION_DESCRIPTION"]
        self.rows = [
            ["BRCA1", "ENST00000357654.9", "COSO1", "missense_variant"],
            ["BRCA1", "ENST00000357654.9", "COSO1", "missense_variant"],
            ["BRCA1", "OTHER", "COSO2", "frameshift_variant"],
            ["TP53", "ENST_OTHER", "COSO1", "missense_variant"],
        ]

    def source(self, rows=None, columns=None, unsafe=False, corrupt=False):
        text = io.StringIO()
        writer = csv.writer(text, delimiter="\t", lineterminator="\n")
        writer.writerow(self.columns if columns is None else columns)
        writer.writerows(self.rows if rows is None else rows)
        payload = b"not gzip" if corrupt else gzip.compress(text.getvalue().encode(), mtime=0)
        with tarfile.open(self.archive, "w") as bundle:
            name = "Cosmic_MutantCensus_v101_GRCh38.tsv.gz"
            member = tarfile.TarInfo(("../" if unsafe else "nested/") + name)
            member.size = len(payload)
            bundle.addfile(member, io.BytesIO(payload))
        self.manifest.write_text(json.dumps({"release": "v101", "genome": "GRCh38", "files": [
            {"product": "mutantcensus", "filename": self.archive.name, "sha256": prep.sha256(self.archive)}]}))

    def run_preparation(self):
        return prep.prepare(self.archive, self.out, self.manifest)

    def test_preserves_duplicate_indel_and_other_transcript_rows_without_tissue_claim(self):
        self.source()
        report = self.run_preparation()
        self.assertEqual(report["counts"], {"source_rows": 4, "brca1_rows": 3})
        self.assertEqual(report["transcripts"], {"ENST00000357654.9": 2, "OTHER": 1})
        self.assertFalse(report["breast_filter_applied"])
        self.assertFalse(report["observed_cosmic_variant_overlay_reproduced"])
        self.assertTrue(report["classification_required_for_breast_selection"])
        with (self.out / "brca1_source_rows.tsv").open() as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        self.assertEqual(rows, [self.columns] + self.rows[:3])

    def test_wrong_hash_does_not_publish(self):
        self.source()
        with self.archive.open("ab") as stream:
            stream.write(b"modified")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            self.run_preparation()
        self.assertFalse(self.out.exists())

    def test_existing_output_is_preserved(self):
        self.source()
        self.run_preparation()
        original = (self.out / "cosmic_brca1_preparation_summary.json").read_bytes()
        with self.assertRaises(FileExistsError):
            self.run_preparation()
        self.assertEqual((self.out / "cosmic_brca1_preparation_summary.json").read_bytes(), original)

    def test_no_brca1_is_an_error(self):
        self.source(rows=[self.rows[-1]])
        with self.assertRaisesRegex(ValueError, "No exact BRCA1"):
            self.run_preparation()
        self.assertFalse(self.out.exists())

    def test_rejects_traversal_without_extracting(self):
        self.source(unsafe=True)
        with self.assertRaisesRegex(ValueError, "unsafe TAR"):
            self.run_preparation()
        self.assertFalse(self.out.exists())
        self.assertFalse((self.root / "Cosmic_MutantCensus_v101_GRCh38.tsv.gz").exists())

    def test_corrupt_gzip_does_not_publish(self):
        self.source(corrupt=True)
        with self.assertRaises(OSError):
            self.run_preparation()
        self.assertFalse(self.out.exists())

    def test_malformed_tsv_does_not_publish(self):
        self.source(rows=[["BRCA1", "only two fields"]])
        with self.assertRaisesRegex(ValueError, "Malformed TSV"):
            self.run_preparation()
        self.assertFalse(self.out.exists())


if __name__ == "__main__":
    unittest.main()
