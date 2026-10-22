"""Independent strand/boundary checks and synthetic tissue-join integration cases."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "01_data_downloading_and_preprocessing"))
import prepare_cosmic_breast_variants as breast


class MappingFixture:
    def setUp(self):
        self.cds = "ATGCCGAAATAA"
        self.mapping = {"strand": "-", "chromosome": "17", "cds_nt": 12, "segments": [
            {"genome_start": 200, "genome_end": 205, "cds_start": 1, "cds_end": 6},
            {"genome_start": 100, "genome_end": 105, "cds_start": 7, "cds_end": 12}]}
        self.positions = breast.position_map(self.mapping)

    def row(self, position=205, reference="T", alternative="C"):
        return {"CHROMOSOME": "17", "GENOME_START": str(position), "GENOME_STOP": str(position),
                "GENOMIC_WT_ALLELE": reference, "GENOMIC_MUT_ALLELE": alternative,
                "MUTATION_CDS": "c.999C>T", "TRANSCRIPT_ACCESSION": "DIFFERENT_SOURCE_TRANSCRIPT"}


class MappingTests(MappingFixture, unittest.TestCase):
    def test_minus_strand_splice_boundaries_and_terminal_stop(self):
        self.assertEqual({p: self.positions[p] for p in [205, 200, 105, 100]}, {205: 1, 200: 6, 105: 7, 100: 12})
        for genomic, expected in [(205, 1), (105, 7), (100, 12)]:
            result = breast.map_variant(self.row(genomic), self.mapping, self.positions, self.cds)
            self.assertEqual((result["status"], result["cds_position"], result["cds_ref"], result["cds_alt"]),
                             ("mapped", expected, "A", "G"))
        self.assertEqual(breast.map_variant(self.row(150), self.mapping, self.positions, self.cds)["status"], "outside_target_cds")

    def test_plus_strand_does_not_complement(self):
        mapping = {"strand": "+", "chromosome": "17", "cds_nt": 3, "segments": [
            {"genome_start": 100, "genome_end": 102, "cds_start": 1, "cds_end": 3}]}
        result = breast.map_variant(self.row(100, "A", "C"), mapping, breast.position_map(mapping), "ATG")
        self.assertEqual(result["cds_substitution"], "1A>C")

    def test_mismatch_and_hgvs_conflict_are_not_accepted(self):
        self.assertEqual(breast.map_variant(self.row(205, "C", "A"), self.mapping, self.positions, self.cds)["status"], "reference_mismatch")
        row = self.row()
        row["HGVSG"] = "17:g.204T>C"
        self.assertEqual(breast.map_variant(row, self.mapping, self.positions, self.cds)["status"], "genomic_annotation_conflict")

    def test_non_snv_is_not_mapped(self):
        self.assertEqual(breast.map_variant(self.row(205, "TT", "C"), self.mapping, self.positions, self.cds)["status"], "not_snv_or_unknown_allele")

    def test_real_frozen_map_covers_5592_positions_and_terminal_stop(self):
        mapping = json.loads(breast.DEFAULT_MAP.read_text())
        positions = breast.position_map(mapping)
        self.assertEqual(len(positions), 5592)
        self.assertEqual(positions[43124096], 1)
        self.assertEqual(positions[43045678], 5592)
        self.assertEqual(mapping["source_gtf"]["official_md5"], "14a867b82917c8c3006838c3a5053a3e")


class IntegrationTests(MappingFixture, unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = self.root / "processed"
        self.fasta = self.root / "input.fa"
        sequence = "G" + self.cds + "C"
        self.fasta.write_text(">TARGET|BRCA1\n" + sequence + "\n")
        self.mapping.update(transcript="TARGET", cds_sha256=hashlib.sha256(self.cds.encode()).hexdigest(),
                            cdsplus_sha256=hashlib.sha256(sequence.encode()).hexdigest(), source_gtf={"sha256": "synthetic"})
        self.mapfile = self.root / "map.json"
        self.mapfile.write_text(json.dumps(self.mapping))
        base = {**self.row(), "GENE_SYMBOL": "BRCA1", "COSMIC_PHENOTYPE_ID": "COSO1",
                "COSMIC_SAMPLE_ID": "COSS1", "GENOMIC_MUTATION_ID": "COSV1"}
        self.observations = [base.copy(), base.copy(), {**base, "GENOME_START": "150", "GENOME_STOP": "150"},
                             {**base, "GENOMIC_WT_ALLELE": "TT"}, {**base, "COSMIC_PHENOTYPE_ID": "COSO2"},
                             {**base, "GENE_SYMBOL": "TP53"}]
        self.classification = [{"COSMIC_PHENOTYPE_ID": "COSO1", "PRIMARY_SITE": "breast"},
                               {"COSMIC_PHENOTYPE_ID": "COSO2", "PRIMARY_SITE": "lung"}]

    def write_sources(self):
        records = []
        for product, name, rows in [("mutantcensus", "MutantCensus", self.observations),
                                     ("classification", "Classification", self.classification)]:
            text = io.StringIO()
            writer = csv.DictWriter(text, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            payload = gzip.compress(text.getvalue().encode(), mtime=0)
            path = self.root / ("Cosmic_" + name + "_Tsv_v101_GRCh38.tar")
            with tarfile.open(path, "w") as bundle:
                member = tarfile.TarInfo("Cosmic_" + name + "_v101_GRCh38.tsv.gz")
                member.size = len(payload)
                bundle.addfile(member, io.BytesIO(payload))
            records.append({"product": product, "filename": path.name, "sha256": breast.sha256(path)})
        self.manifest = self.root / "sources.json"
        self.manifest.write_text(json.dumps({"release": "v101", "genome": "GRCh38", "files": records}))

    def run_preparation(self):
        self.write_sources()
        return breast.prepare(self.root, self.fasta, self.out, self.manifest, self.mapfile, breast.sha256(self.fasta))

    def test_join_mapping_and_duplicate_observation_accounting(self):
        report = self.run_preparation()
        self.assertTrue(report["preparation_passed"])
        self.assertEqual(report["counts"]["source_rows"], 6)
        self.assertEqual(report["counts"]["brca1_rows"], 5)
        self.assertEqual(report["counts"]["breast_brca1_rows"], 4)
        self.assertEqual(report["counts"]["mapped"], 2)
        self.assertEqual(report["unique_mapped_snvs"], 1)
        with (self.out / "brca1_breast_target_snvs.tsv").open() as stream:
            variants = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(variants[0]["cds_substitution"], "1A>G")
        self.assertEqual((variants[0]["source_observations"], variants[0]["unique_samples"]), ("2", "1"))
        self.assertFalse(report["observed_cosmic_variant_overlay_reproduced"])

    def test_unknown_phenotype_fails_with_audit_output(self):
        self.observations[0]["COSMIC_PHENOTYPE_ID"] = "UNKNOWN"
        report = self.run_preparation()
        self.assertFalse(report["preparation_passed"])
        self.assertIn("unresolved_phenotype", report["attention_reasons"])
        self.assertTrue((self.out / "brca1_observation_audit.tsv").exists())

    def test_reference_mismatch_fails_despite_other_valid_variants(self):
        self.observations[0]["GENOMIC_WT_ALLELE"] = "C"
        self.observations[0]["GENOMIC_MUT_ALLELE"] = "A"
        report = self.run_preparation()
        self.assertFalse(report["preparation_passed"])
        self.assertIn("reference_mismatch", report["attention_reasons"])

    def test_conflicting_classification_rejected(self):
        self.classification.append({"COSMIC_PHENOTYPE_ID": "COSO1", "PRIMARY_SITE": "lung"})
        with self.assertRaisesRegex(ValueError, "Conflicting tissue"):
            self.run_preparation()
        self.assertFalse(self.out.exists())

    def test_wrong_reference_cannot_publish_output(self):
        self.fasta.write_text(self.fasta.read_text().replace(self.cds, "C" + self.cds[1:]))
        with self.assertRaisesRegex(ValueError, "reference sequence"):
            self.run_preparation()
        self.assertFalse(self.out.exists())

    def test_existing_output_not_overwritten(self):
        self.run_preparation()
        with self.assertRaises(FileExistsError):
            self.run_preparation()


if __name__ == "__main__":
    unittest.main()
