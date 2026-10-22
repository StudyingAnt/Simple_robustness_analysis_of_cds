#!/usr/bin/env python3
"""Compare audited BRCA1 graph states with mapped breast-cancer SNVs without resimulation."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time

from validate_manuscript_results import digest, require, save_plot
from run_supplementary_cases import audit_graph
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "01_data_downloading_and_preprocessing"
TARGET = "ENST00000357654.9"
GRAPH_FILES = ("s12_graph_nodes.csv", "s12_graph_edges.csv", "trajectory_counts.csv", "s12_reference.cdsplus.fa")
VARIANT_FILE = "brca1_breast_target_snvs.tsv"


def verified_files(folder, report, names):
    entries = report.get("outputs", [])
    require(len({r["filename"] for r in entries}) == len(entries), "Duplicate output manifest entries")
    by_name = {r["filename"]: r["sha256"] for r in entries}
    result = {}
    for name in names:
        require(name in by_name and Path(name).name == name, "Missing artifact hash: " + name)
        path = folder / name
        require(path.resolve().parent == folder.resolve() and digest(path) == by_name[name], "Artifact changed: " + str(path))
        result[name] = by_name[name]
    return result


def check_graph_manifest(folder):
    report = json.loads((folder / "supplementary_summary.json").read_text())
    require(report.get("stage") == "s12" and report.get("simulation_passed") is True
            and report.get("graph_audit_passed") is True, "Graph report has not passed")
    require(report.get("selected_transcripts") == [TARGET] and report.get("n_trajectories") == 1000
            and report.get("complete_paths") == 1000 and report.get("incomplete_paths") == 0
            and report.get("bounded_software_check") is False, "Require the complete 1,000-path BRCA1 graph")
    require(report.get("signatures") == ["SBS2"] and report.get("seed") == 42, "Unexpected graph model/seed")
    require(report.get("signature_version") == "3.3.1" and report.get("signature_reference") == "GRCh37"
            and report.get("signature_sha256") == "9d9c212cc14aa660da1b7fd9e60e696c8e2c885a7606c27d20ad8997954554ea",
            "Graph signature profile differs from the pinned simulation")
    criterion = report.get("criterion", {})
    require(all(criterion.get(k) == v for k, v in dict(matrix="BLOSUM80", threshold=20,
            mode="global", gap_open=-10, gap_extend=-.5, stop_policy="score").items()), "Graph criterion differs")
    mapping = json.loads((DATA / "brca1_gencode_v40_coordinates.json").read_text())
    require(report.get("reference_sequence_sha256") == mapping["cdsplus_sha256"], "Unexpected graph reference")
    return report, verified_files(folder, report, GRAPH_FILES)


def check_variant_manifest(folder):
    report = json.loads((folder / "cosmic_breast_variants_summary.json").read_text())
    require(report.get("stage") == "cosmic_brca1_breast_coordinate_mapping" and report.get("preparation_passed") is True
            and report.get("attention_reasons") == [], "Variant preparation has not passed")
    require(report.get("release") == "v101" and report.get("genome") == "GRCh38"
            and report.get("target_transcript") == TARGET and report.get("target_cds_nt") == 5592
            and report.get("breast_filter_applied") is True and report.get("snv_filter_applied") is True,
            "Unexpected variant selection/target")
    manifest = json.loads((DATA / "cosmic_v101_grch38.json").read_text())
    inputs = {r["product"]: r["sha256"] for r in report["inputs"]}
    require(inputs == {r["product"]: r["sha256"] for r in manifest["files"]}, "Variant source archive fingerprints differ")
    require(report["source_manifest_sha256"] == digest(DATA / "cosmic_v101_grch38.json")
            and report["coordinate_map_sha256"] == digest(DATA / "brca1_gencode_v40_coordinates.json"),
            "Variant source/coordinate-map manifest differs")
    return report, verified_files(folder, report, (VARIANT_FILE, "brca1_observation_audit.tsv"))


def discover(root, pattern, checker):
    """Select only successful verified inputs; equivalent duplicate runs are harmless."""
    valid = []
    rejected = []
    for folder in sorted(Path(root).glob(pattern)):
        if not folder.is_dir():
            continue
        try:
            _, fingerprints = checker(folder)
            valid.append((folder, fingerprints))
        except (ValueError, KeyError, OSError, TypeError) as error:
            rejected.append(str(folder) + ": " + str(error))
    if not valid:
        raise ValueError("No successful verified input under " + str(root) + "; " + "; ".join(rejected))
    fingerprints = {json.dumps(hashes, sort_keys=True) for _, hashes in valid}
    require(len(fingerprints) == 1, "Multiple scientifically different inputs found; provide an explicit input directory")
    return valid[-1][0]


def read_reference(path):
    lines = path.read_text().splitlines()
    require(lines and lines[0].startswith(">" + TARGET + "|")
            and not any(line.startswith(">") for line in lines[1:]), "Unexpected graph reference FASTA")
    return "".join(line.strip() for line in lines[1:]).upper()


def read_graph(folder, report):
    reference = read_reference(folder / "s12_reference.cdsplus.fa")
    require(hashlib.sha256(reference.encode()).hexdigest() == report["reference_sequence_sha256"], "Reference sequence hash differs")
    with (folder / "s12_graph_nodes.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    nodes = {}
    for row in rows:
        require(row["node"] not in nodes and row["transcript_id"] == TARGET, "Duplicate/wrong-transcript graph node")
        for key in ("hamming_distance",):
            row[key] = int(row[key])
        for key in ("blosum_distance", "x", "y"):
            row[key] = float(row[key])
            require(math.isfinite(row[key]), "Non-finite graph value")
        require(abs(math.hypot(row["x"], row["y"]) - row["hamming_distance"]) < 1e-7, "Graph radial layout is inconsistent")
        nodes[row["node"]] = row
    with (folder / "s12_graph_edges.csv").open(newline="") as stream:
        edges = list(csv.DictReader(stream))
    for row in edges:
        row["path"], row["step"] = int(row["path"]), int(row["step"])
    with (folder / "trajectory_counts.csv").open(newline="") as stream:
        paths = list(csv.DictReader(stream))
    for row in paths:
        row["path"], row["blosum_neutral_steps"] = int(row["path"]), int(row["blosum_neutral_steps"])
        row["terminal_distance"] = float(row["terminal_distance"])
    require(len(nodes) == report["nodes"] and len(edges) == report["edges"] and len(paths) == 1000, "Graph counts differ")
    audit = audit_graph(paths, nodes, edges, reference)
    return reference, nodes, edges, audit


def read_variants(path, reference, report):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    by_change = {}
    for row in rows:
        position = int(row["cds_position"])
        old, new = row["cds_ref"], row["cds_alt"]
        change = str(position) + old + ">" + new
        require(0 < position < len(reference) - 1 and reference[position] == old
                and old in "ACGT" and new in "ACGT" and len(old) == len(new) == 1 and old != new,
                "Variant coordinate/reference is inconsistent")
        require(row["target_transcript"] == TARGET and row["cds_substitution"] == change
                and change not in by_change, "Wrong transcript, substitution label or duplicate variant")
        observations, samples = int(row["source_observations"]), int(row["unique_samples"])
        require(observations > 0 and 0 <= samples <= observations, "Invalid observation/sample counts")
        by_change[change] = row
    require(len(rows) == report["unique_mapped_snvs"]
            and sum(int(r["source_observations"]) for r in rows) == report["counts"]["mapped"],
            "Prepared variant counts differ")
    return by_change


def annotate(nodes, observed):
    annotated = []
    per_variant = {change: dict(exact_single_snv_nodes=0, component_match_nodes=0) for change in observed}
    counts = Counter()
    for key in sorted(nodes):
        node = nodes[key]
        changes = set(filter(None, node["cds_substitutions"].split(";")))
        require(len(changes) == node["hamming_distance"], "Invalid graph substitutions")
        matched = sorted(changes.intersection(observed))
        exact = len(changes) == 1 and len(matched) == 1
        multi = len(changes) > 1 and bool(matched)
        counts["nodes"] += 1
        counts["exact_single_snv_nodes"] += exact
        counts["multi_mutation_nodes_with_reported_component"] += multi
        counts["nodes_with_any_reported_component"] += bool(matched)
        for change in matched:
            per_variant[change]["component_match_nodes"] += 1
            per_variant[change]["exact_single_snv_nodes"] += exact
        status = ("reference" if not changes else "single_snv_allele_match" if exact
                  else "component_allele_match_only" if matched else "no_reported_component_match")
        annotated.append(dict(node, simulation_observed_cosmic_status=node["observed_cosmic_status"],
            observed_cosmic_status=status, cosmic_exact_single_snv_match=exact,
            cosmic_matched_component_count=len(matched), cosmic_matched_components=";".join(matched),
            cosmic_all_components_reported=bool(changes) and len(matched) == len(changes),
            observed_genotype_confirmed=False))
    counts["reported_variants_matching_any_component"] = sum(r["component_match_nodes"] > 0 for r in per_variant.values())
    counts["reported_variants_matching_single_snv_nodes"] = sum(r["exact_single_snv_nodes"] > 0 for r in per_variant.values())
    return annotated, per_variant, dict(counts)


def write_table(path, rows, fields, delimiter="\t"):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot(stage, rows, edges, software_check=False):
    nodes = {r["node"]: r for r in rows}
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.add_collection(LineCollection([[(nodes[e["source"]]["x"], nodes[e["source"]]["y"]),
                                     (nodes[e["target"]]["x"], nodes[e["target"]]["y"])] for e in edges],
                                    color="0.7", alpha=.18, linewidth=.35, rasterized=True))
    scatter = ax.scatter([r["x"] for r in rows], [r["y"] for r in rows],
                         c=[r["blosum_distance"] for r in rows], cmap="viridis", s=9, linewidths=0,
                         rasterized=True, zorder=2)
    fig.colorbar(scatter, ax=ax, label="BLOSUM80 alignment-score distance", shrink=.8)
    matched = [r for r in rows if r["cosmic_exact_single_snv_match"]]
    ax.scatter([r["x"] for r in matched], [r["y"] for r in matched], color="#d62728", edgecolors="white",
               linewidths=.45, s=38, zorder=4)
    ax.scatter([0], [0], marker="*", color="black", s=65, zorder=5)
    ax.legend(handles=[Line2D([], [], marker="o", linestyle="", color="#d62728",
                            label=f"Single-SNV allele match (n={len(matched)})"),
                       Line2D([], [], marker="*", linestyle="", color="black", label="Reference")],
              loc="upper right", frameon=False, fontsize=8)
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set(xlabel="Radius = reference DNA Hamming distance", ylabel="Deterministic angular layout",
           title="BRCA1-203: SBS2 simulation and COSMIC breast SNV matches")
    if software_check:
        ax.set_title("SOFTWARE CHECK: synthetic alleles on the BRCA1 graph")
    fig.text(.5, .015, "Red marks allele agreement in a one-substitution state; complete tumour genotypes are not established.",
             ha="center", fontsize=7)
    fig.subplots_adjust(bottom=.14, top=.9)
    save_plot(fig, stage, "figureS12_cosmic_single_snv_matches")


def run(graph_dir, variants_dir, output):
    started = time.perf_counter()
    graph_dir, variants_dir, output = Path(graph_dir).resolve(), Path(variants_dir).resolve(), Path(output).resolve()
    require(not output.exists(), "Choose a new annotation output directory")
    for source in (graph_dir, variants_dir):
        require(output != source and source not in output.parents, "Keep annotation outputs outside source directories")
    graph_report, graph_hashes = check_graph_manifest(graph_dir)
    variants_report, variant_hashes = check_variant_manifest(variants_dir)
    reference, nodes, edges, audit = read_graph(graph_dir, graph_report)
    require(graph_report["input_sha256"] == variants_report["human_fasta_sha256"], "Human reference inputs differ")
    observed = read_variants(variants_dir / VARIANT_FILE, reference, variants_report)
    rows, per_variant, counts = annotate(nodes, observed)
    output.parent.mkdir(parents=True, exist_ok=True)
    source_report_hashes = {str(graph_dir / "supplementary_summary.json"): digest(graph_dir / "supplementary_summary.json"),
                            str(variants_dir / "cosmic_breast_variants_summary.json"): digest(variants_dir / "cosmic_breast_variants_summary.json")}
    with tempfile.TemporaryDirectory(prefix=".brca1-observed-", dir=output.parent) as temp:
        stage = Path(temp) / "results"
        stage.mkdir()
        write_table(stage / "s12_annotated_nodes.tsv", rows, list(rows[0]))
        comparisons = [dict(observed[key], **per_variant[key]) for key in sorted(observed, key=lambda s: (int(re.match(r'\d+', s)[0]), s))]
        require(comparisons, "Observed variant set is empty")
        write_table(stage / "s12_observed_variant_matches.tsv", comparisons, list(comparisons[0]))
        software_check = variants_report.get("bounded_software_check", False) is True
        plot(stage, rows, edges, software_check=software_check)
        caption = ("BRCA1-203 (ENST00000357654.9), 1,000 SBS2 trajectories, BLOSUM80 distance threshold 20. "
                   "Nodes are simulated DNA states and edges are chronological single-substitution events. "
                   "Node colour shows BLOSUM80 alignment-score distance; radius shows reference DNA Hamming distance. "
                   "Red points are states with exactly one reference-relative CDS substitution matching a mapped "
                   "COSMIC v101 GRCh38 breast-cancer SNV. Matches indicate reported allele agreement and do not establish "
                   "complete tumour genotypes, causation by APOBEC, or driver status. Multi-substitution component matches "
                   "are recorded separately and are not coloured as observed genotypes. The historical release and "
                   "observation-to-node rule remain unconfirmed; this is a documented reconstruction.\n")
        if software_check:
            caption = "SOFTWARE CHECK WITH SYNTHETIC ALLELES; NOT A COSMIC OBSERVATION RESULT.\n" + caption
        (stage / "figureS12_caption.txt").write_text(caption)
        report = dict(stage="s12_observed_snv_annotation", annotation_passed=True, fresh_simulation=False,
            original_figure_reproduced=False, historical_observation_rule_recovered=False,
            regenerated_allele_match_panel_created=True,
            bounded_software_check=software_check,
            observed_snv_matching_performed=True, observed_genotypes_confirmed=False,
            observed_cosmic_variant_overlay_reproduced=False,
            target_transcript=TARGET, n_trajectories=1000, cosmic_release="v101", cosmic_genome="GRCh38",
            distinct_observed_snvs=len(observed), source_observations=variants_report["counts"]["mapped"],
            matching_rule="Exactly one reference-relative substitution and that allele appears in the mapped breast SNV set.",
            counts=counts, graph_audit=audit, graph_directory=str(graph_dir), variants_directory=str(variants_dir),
            graph_file_hashes=graph_hashes, variant_file_hashes=variant_hashes,
            input_report_hashes=source_report_hashes, script_sha256=digest(__file__),
            interpretation="A regenerated figure with explicit allele matching. Matching individual variants does not prove that a simulated multi-mutation genotype was observed. Historical artwork/annotation semantics are not claimed reproduced.",
            elapsed_seconds=time.perf_counter()-started,
            outputs=[dict(filename=p.name, sha256=digest(p)) for p in sorted(stage.iterdir())])
        require(verified_files(graph_dir, graph_report, GRAPH_FILES) == graph_hashes
                and verified_files(variants_dir, variants_report, (VARIANT_FILE, "brca1_observation_audit.tsv")) == variant_hashes
                and all(digest(Path(path)) == sha for path, sha in source_report_hashes.items()), "Inputs changed during annotation")
        (stage / "s12_observed_annotation_summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        require(not output.exists(), "Output appeared during annotation")
        stage.rename(output)
    print(json.dumps({k: report[k] for k in ("stage", "annotation_passed", "fresh_simulation", "distinct_observed_snvs", "counts")}, indent=2))
    print("Report:", output / "s12_observed_annotation_summary.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    graph = parser.add_mutually_exclusive_group(required=True)
    graph.add_argument("--graph-dir", type=Path)
    graph.add_argument("--graph-root", type=Path, help="Discover verified runs/brca1_sbs2_blosum20_n1000_* input")
    variants = parser.add_mutually_exclusive_group(required=True)
    variants.add_argument("--variants-dir", type=Path)
    variants.add_argument("--variants-root", type=Path, help="Discover verified processed/brca1_breast_* input")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    graph_dir = args.graph_dir or discover(args.graph_root, "brca1_sbs2_blosum20_n1000_*", check_graph_manifest)
    variants_dir = args.variants_dir or discover(args.variants_root, "brca1_breast_*", check_variant_manifest)
    run(graph_dir, variants_dir, args.out_dir)


if __name__ == "__main__":
    main()
