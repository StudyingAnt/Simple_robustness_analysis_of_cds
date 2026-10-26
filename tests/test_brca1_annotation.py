"""Allele-match semantics and input integrity without restricted COSMIC observations."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "07_robustness_analysis"))
import annotate_brca1_graph as annotation


class BrCaAnnotationTests(unittest.TestCase):
    def nodes(self):
        result = {}
        for key, changes, distance in [("root", "", 0), ("one", "1A>G", 1),
                                       ("two", "1A>G;2C>T", 2), ("other", "2C>A", 1)]:
            result[key] = dict(node=key, transcript_id=annotation.TARGET, cds_substitutions=changes,
                               hamming_distance=distance, blosum_distance=float(distance),
                               observed_cosmic_status="not_available", x=float(distance), y=0.)
        return result

    def test_component_matches_never_confirm_multi_mutation_genotypes(self):
        rows, variants, counts = annotation.annotate(self.nodes(), {"1A>G": {}, "2C>T": {}})
        by_key = {row["node"]: row for row in rows}
        self.assertTrue(by_key["one"]["cosmic_exact_single_snv_match"])
        self.assertFalse(by_key["two"]["cosmic_exact_single_snv_match"])
        self.assertTrue(by_key["two"]["cosmic_all_components_reported"])
        self.assertEqual(by_key["two"]["observed_cosmic_status"], "component_allele_match_only")
        self.assertTrue(all(row["observed_genotype_confirmed"] is False for row in rows))
        self.assertEqual(counts["exact_single_snv_nodes"], 1)
        self.assertEqual(counts["multi_mutation_nodes_with_reported_component"], 1)
        self.assertEqual(variants["2C>T"]["exact_single_snv_nodes"], 0)
        self.assertEqual(variants["2C>T"]["component_match_nodes"], 1)

    def test_position_match_with_different_alternate_is_not_an_allele_match(self):
        rows, _, counts = annotation.annotate(self.nodes(), {"2C>T": {}})
        self.assertFalse(next(r for r in rows if r["node"] == "other")["cosmic_exact_single_snv_match"])
        self.assertEqual(counts["exact_single_snv_nodes"], 0)

    def variant_table(self, root, change="1A>G", observations=2):
        path = root / "variants.tsv"
        row = dict(target_transcript=annotation.TARGET, cds_position=1, cds_ref="A", cds_alt="G",
                   cds_substitution=change, source_observations=observations, unique_samples=1)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        return path

    def test_variants_checked_against_reference_and_prepared_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.variant_table(Path(directory))
            report = dict(unique_mapped_snvs=1, counts=dict(mapped=2))
            self.assertEqual(set(annotation.read_variants(path, "GACG", report)), {"1A>G"})
            with self.assertRaisesRegex(ValueError, "coordinate/reference"):
                annotation.read_variants(path, "GCCG", report)
            report["counts"]["mapped"] = 3
            with self.assertRaisesRegex(ValueError, "counts differ"):
                annotation.read_variants(path, "GACG", report)

    def test_variant_label_is_not_trusted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.variant_table(Path(directory), change="2C>T")
            with self.assertRaisesRegex(ValueError, "substitution label"):
                annotation.read_variants(path, "GACG", dict(unique_mapped_snvs=1, counts=dict(mapped=2)))

    def test_modified_artifact_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "input.csv"
            path.write_text("original")
            report = dict(outputs=[dict(filename=path.name, sha256=annotation.digest(path))])
            self.assertEqual(annotation.verified_files(root, report, (path.name,))[path.name], annotation.digest(path))
            path.write_text("changed")
            with self.assertRaisesRegex(ValueError, "Artifact changed"):
                annotation.verified_files(root, report, (path.name,))

    def test_discovery_requires_explicit_selection_when_scientific_inputs_differ(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ["run_a", "run_b"]:
                (root / name).mkdir()
            with self.assertRaisesRegex(ValueError, "scientifically different"):
                annotation.discover(root, "run_*", lambda path: ({}, {"data": path.name}))

    def test_discovery_reuses_only_equivalent_verified_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ["run_a", "run_b"]:
                (root / name).mkdir()
            chosen = annotation.discover(root, "run_*", lambda path: ({}, {"data": "same hash"}))
            self.assertEqual(chosen.name, "run_b")

    def test_failed_preparation_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "cosmic_breast_variants_summary.json").write_text(json.dumps(
                dict(stage="cosmic_brca1_breast_coordinate_mapping", preparation_passed=False)))
            with self.assertRaisesRegex(ValueError, "has not passed"):
                annotation.check_variant_manifest(root)

    def test_plot_accepts_zero_real_matches_without_inventing_highlights(self):
        with tempfile.TemporaryDirectory() as directory:
            rows, _, counts = annotation.annotate(self.nodes(), {"3G>A": {}})
            self.assertEqual(counts["exact_single_snv_nodes"], 0)
            annotation.plot(Path(directory), rows, [])
            self.assertTrue((Path(directory) / "figureS12_cosmic_single_snv_matches.png").exists())

    def test_existing_output_is_preserved_before_input_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, "new annotation output"):
                annotation.run(root / "missing_graph", root / "missing_variants", root)


if __name__ == "__main__":
    unittest.main()
