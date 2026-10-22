#!/usr/bin/env python3
"""Select breast observations and map genomic BRCA1 SNVs to the simulated CDS."""
import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile

from build_brca1_coordinate_map import DEFAULT_MAP
from download_cosmic_variants import MANIFEST, sha256

HUMAN_FASTA_SHA256 = "c2f30678cb9f2252bd4717d1f0b8fef9e9575b4be843aa64478c4b1e0422a751"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@contextmanager
def table(archive, record):
    name = record["filename"].replace("_Tsv_", "_").removesuffix(".tar") + ".tsv.gz"
    with tarfile.open(archive, "r:*") as bundle:
        members = bundle.getmembers()
        if any(PurePosixPath(m.name).is_absolute() or ".." in PurePosixPath(m.name).parts
               or m.issym() or m.islnk() for m in members):
            raise ValueError("Unsafe archive member.")
        matches = [m for m in members if m.isfile() and PurePosixPath(m.name).name == name]
        if len(matches) != 1:
            raise ValueError("Expected exactly one table: " + name)
        with bundle.extractfile(matches[0]) as raw, gzip.GzipFile(fileobj=raw) as decoded:
            with io.TextIOWrapper(decoded, encoding="utf-8-sig", newline="") as text:
                reader = csv.DictReader(text, delimiter="\t")
                fields = reader.fieldnames or []
                if not fields or len(fields) != len(set(fields)):
                    raise ValueError("Missing/duplicate table columns.")
                yield reader


def check_row(row):
    if None in row or any(value is None for value in row.values()):
        raise ValueError("Malformed TSV row.")


def target_sequence(fasta, mapping):
    records = []
    selected = False
    chunks = []
    with Path(fasta).open() as stream:
        for line in stream:
            if line.startswith(">"):
                if selected:
                    records.append("".join(chunks))
                selected = line[1:].strip().split("|", 1)[0].split()[0] == mapping["transcript"]
                chunks = []
            elif selected:
                chunks.append(line.strip().upper())
        if selected:
            records.append("".join(chunks))
    if len(records) != 1:
        raise ValueError("Target transcript missing or duplicated in CDS-plus FASTA.")
    sequence = records[0]
    cds = sequence[1:-1]
    if (len(cds) != mapping["cds_nt"] or hashlib.sha256(cds.encode()).hexdigest() != mapping["cds_sha256"]
            or hashlib.sha256(sequence.encode()).hexdigest() != mapping["cdsplus_sha256"]):
        raise ValueError("Target reference sequence does not match the simulated CDS-plus sequence.")
    return cds


def position_map(mapping):
    if mapping["strand"] not in ("+", "-"):
        raise ValueError("Invalid target strand.")
    positions = {}
    for segment in mapping["segments"]:
        start, end = segment["genome_start"], segment["genome_end"]
        if start < 1 or end < start or segment["cds_end"] - segment["cds_start"] != end - start:
            raise ValueError("Invalid coding interval.")
        ordered = range(start, end + 1) if mapping["strand"] == "+" else range(end, start - 1, -1)
        for offset, genomic in enumerate(ordered):
            if genomic in positions:
                raise ValueError("Overlapping coding intervals.")
            positions[genomic] = segment["cds_start"] + offset
    if sorted(positions.values()) != list(range(1, mapping["cds_nt"] + 1)):
        raise ValueError("Non-contiguous target CDS coordinates.")
    return positions


def map_variant(row, mapping, positions, cds):
    chrom = row["CHROMOSOME"].strip().removeprefix("chr")
    if chrom != mapping["chromosome"]:
        return {"status": "other_chromosome"}
    reference, alternative = row["GENOMIC_WT_ALLELE"].strip().upper(), row["GENOMIC_MUT_ALLELE"].strip().upper()
    if len(reference) != 1 or len(alternative) != 1 or reference not in "ACGT" or alternative not in "ACGT" or reference == alternative:
        return {"status": "not_snv_or_unknown_allele"}
    try:
        start, end = int(row["GENOME_START"]), int(row["GENOME_STOP"])
    except ValueError:
        return {"status": "invalid_genomic_coordinates"}
    if start < 1 or start != end:
        return {"status": "not_single_position"}
    genomic_hgvs = row.get("HGVSG", "").strip()
    match = re.fullmatch(r"[^:]+:g\.(\d+)([ACGT])>([ACGT])", genomic_hgvs)
    if match and (int(match[1]) != start or match[2] != reference or match[3] != alternative):
        return {"status": "genomic_annotation_conflict"}
    if start not in positions:
        return {"status": "outside_target_cds"}
    coordinate = positions[start]
    if mapping["strand"] == "-":
        reference, alternative = reference.translate(COMPLEMENT), alternative.translate(COMPLEMENT)
    if cds[coordinate - 1] != reference:
        return {"status": "reference_mismatch", "cds_position": coordinate,
                "cds_ref": reference, "cds_alt": alternative}
    context = cds[max(0, coordinate - 2):coordinate + 1]
    return {"status": "mapped", "cds_position": coordinate, "cds_ref": reference,
            "cds_alt": alternative, "cds_substitution": str(coordinate) + reference + ">" + alternative,
            "target_cds_context": context if len(context) == 3 and coordinate > 1 else "boundary"}


def prepare(raw_dir, fasta, out_dir, manifest_path=MANIFEST, map_path=DEFAULT_MAP,
            expected_fasta_sha256=HUMAN_FASTA_SHA256):
    raw_dir, fasta, out_dir = Path(raw_dir), Path(fasta), Path(out_dir)
    if out_dir.exists():
        raise FileExistsError("Choose a new output directory: " + str(out_dir))
    manifest = json.loads(Path(manifest_path).read_text())
    mapping = json.loads(Path(map_path).read_text())
    records = {r["product"]: r for r in manifest["files"]}
    inputs = []
    for key in ("mutantcensus", "classification"):
        archive = raw_dir / records[key]["filename"]
        digest = sha256(archive)
        if digest != records[key]["sha256"]:
            raise ValueError("Official archive SHA256 mismatch: " + archive.name)
        inputs.append({"product": key, "path": str(archive.resolve()), "sha256": digest})
    if sha256(fasta) != expected_fasta_sha256:
        raise ValueError("Human FASTA SHA256 mismatch.")
    cds = target_sequence(fasta, mapping)
    positions = position_map(mapping)
    classes = {}
    classification_columns = []
    with table(raw_dir / records["classification"]["filename"], records["classification"]) as reader:
        classification_columns = reader.fieldnames
        if not {"COSMIC_PHENOTYPE_ID", "PRIMARY_SITE"}.issubset(reader.fieldnames):
            raise ValueError("Classification lacks phenotype ID/primary site.")
        for row in reader:
            check_row(row)
            key, site = row["COSMIC_PHENOTYPE_ID"].strip(), row["PRIMARY_SITE"].strip().casefold()
            if not key:
                raise ValueError("Empty Classification phenotype ID.")
            # Missing tissue only matters if this phenotype occurs in the BRCA1 observations.
            site = site or None
            if key in classes and classes[key] != site:
                raise ValueError("Conflicting tissue classification for a phenotype ID.")
            classes[key] = site
    if "breast" not in classes.values():
        raise ValueError("Classification has no exact primary-site value 'breast'.")
    counts = Counter()
    sites = Counter()
    transcripts = Counter()
    audit_rows = []
    accepted = []
    with table(raw_dir / records["mutantcensus"]["filename"], records["mutantcensus"]) as reader:
        required = {"GENE_SYMBOL", "COSMIC_PHENOTYPE_ID", "TRANSCRIPT_ACCESSION", "COSMIC_SAMPLE_ID",
                    "CHROMOSOME", "GENOME_START", "GENOME_STOP", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"}
        if not required.issubset(reader.fieldnames):
            raise ValueError("Mutation table is missing required genomic/sample columns.")
        source_fields = list(reader.fieldnames)
        for row in reader:
            check_row(row)
            counts["source_rows"] += 1
            if row["GENE_SYMBOL"].strip() != "BRCA1":
                continue
            counts["brca1_rows"] += 1
            phenotype = row["COSMIC_PHENOTYPE_ID"].strip()
            site = classes.get(phenotype)
            sites[site or "unresolved_phenotype"] += 1
            if site is None:
                result = {"status": "unresolved_phenotype"}
            elif site != "breast":
                result = {"status": "other_primary_site"}
            else:
                counts["breast_brca1_rows"] += 1
                transcripts[row["TRANSCRIPT_ACCESSION"].strip()] += 1
                result = map_variant(row, mapping, positions, cds)
            counts[result["status"]] += 1
            audit_rows.append({**row, "SELECTED_PRIMARY_SITE": site or "", "MAPPING_STATUS": result["status"],
                               "TARGET_CDS_POSITION": result.get("cds_position", ""),
                               "TARGET_REF": result.get("cds_ref", ""), "TARGET_ALT": result.get("cds_alt", "")})
            if result["status"] == "mapped":
                accepted.append((row, result))
    # Check that every BRCA1 observation is assigned exactly one audit status.
    if counts["brca1_rows"] != len(audit_rows) or sum(sites.values()) != counts["brca1_rows"]:
        raise ValueError("Observation accounting failed.")
    failures = []
    for status in ("unresolved_phenotype", "reference_mismatch", "genomic_annotation_conflict", "invalid_genomic_coordinates"):
        if counts[status]:
            failures.append(status)
    if not accepted:
        failures.append("no_mapped_breast_snvs")
    groups = {}
    for row, mapped in accepted:
        key = (mapped["cds_position"], mapped["cds_ref"], mapped["cds_alt"])
        group = groups.setdefault(key, {"mapped": mapped, "rows": 0, "samples": set(), "cosv": set()})
        group["rows"] += 1
        if row["COSMIC_SAMPLE_ID"].strip():
            group["samples"].add(row["COSMIC_SAMPLE_ID"].strip())
        if row.get("GENOMIC_MUTATION_ID", "").strip():
            group["cosv"].add(row["GENOMIC_MUTATION_ID"].strip())
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".cosmic-breast-", dir=out_dir.parent))
    try:
        extras = ["SELECTED_PRIMARY_SITE", "MAPPING_STATUS", "TARGET_CDS_POSITION", "TARGET_REF", "TARGET_ALT"]
        if set(extras).intersection(source_fields):
            raise ValueError("Source columns collide with audit fields.")
        with (staging / "brca1_observation_audit.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=source_fields + extras, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(audit_rows)
        fields = ["target_transcript", "cds_position", "cds_ref", "cds_alt", "cds_substitution", "target_cds_context",
                  "source_observations", "unique_samples", "cosmic_genomic_ids"]
        with (staging / "brca1_breast_target_snvs.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for key, group in sorted(groups.items()):
                variant = {k: group["mapped"][k] for k in fields if k in group["mapped"]}
                variant.update(target_transcript=mapping["transcript"], source_observations=group["rows"],
                               unique_samples=len(group["samples"]), cosmic_genomic_ids=";".join(sorted(group["cosv"])))
                writer.writerow(variant)
        report = {"stage": "cosmic_brca1_breast_coordinate_mapping", "preparation_passed": not failures,
                  "attention_reasons": failures, "release": manifest["release"], "genome": manifest["genome"],
                  "historical_release_confirmed": False, "breast_filter": "Classification PRIMARY_SITE case-insensitive exact 'breast'",
                  "breast_filter_applied": True, "snv_filter_applied": True,
                  "target_transcript": mapping["transcript"], "target_cds_nt": len(cds),
                  "reference_allele_policy": "Genomic forward-strand alleles, complemented on minus-strand target; mismatch excluded and fails preparation.",
                  "counts": dict(counts), "brca1_primary_sites": dict(sorted(sites.items())),
                  "breast_source_transcripts": dict(sorted(transcripts.items())),
                  "unique_mapped_snvs": len(groups), "classification_columns": classification_columns,
                  "inputs": inputs, "human_fasta_sha256": expected_fasta_sha256,
                  "source_manifest_sha256": sha256(manifest_path), "coordinate_map_sha256": sha256(map_path),
                  "coordinate_gtf_sha256": mapping["source_gtf"]["sha256"],
                  "script_sha256": sha256(__file__), "verified_at_utc": datetime.now(timezone.utc).isoformat(),
                  "observed_cosmic_variant_overlay_reproduced": False,
                  "interpretation": "Breast source observations mapped independently of source transcript CDS numbering. Passing preparation establishes a checked SNV set, not observed multi-mutation genotypes, driver status, APOBEC causation, or reproduction of the historical graph annotation rule.",
                  "outputs": [{"filename": name, "sha256": sha256(staging / name)} for name in [
                      "brca1_observation_audit.tsv", "brca1_breast_target_snvs.tsv"]]}
        # Ensure inputs remained unchanged throughout the joins and mapping.
        if any(sha256(item["path"]) != item["sha256"] for item in inputs) or sha256(fasta) != expected_fasta_sha256:
            raise ValueError("An input changed during preparation.")
        (staging / "cosmic_breast_variants_summary.json").write_text(json.dumps(report, indent=2) + "\n")
        if out_dir.exists():
            raise FileExistsError("Output appeared during processing.")
        os.rename(staging, out_dir)
        return report
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True, help="Pinned human CDS-plus FASTA")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.raw_dir, args.input, args.out_dir)
    print(json.dumps({key: report[key] for key in ["stage", "preparation_passed", "attention_reasons",
          "counts", "breast_source_transcripts", "unique_mapped_snvs", "target_transcript"]}, indent=2))
    print("Report:", args.out_dir / "cosmic_breast_variants_summary.json")
    if not report["preparation_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
