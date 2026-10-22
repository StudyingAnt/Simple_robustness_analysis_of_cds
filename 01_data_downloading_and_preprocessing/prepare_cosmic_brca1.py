#!/usr/bin/env python3
"""Verify the frozen COSMIC archive and retain every BRCA1 source row for inspection."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import io
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile

from download_cosmic_variants import MANIFEST, sha256


def prepare(archive, out_dir, manifest_path=MANIFEST):
    archive, out_dir = Path(archive), Path(out_dir)
    if out_dir.exists():
        raise FileExistsError("Choose a new output directory: " + str(out_dir))
    manifest = json.loads(Path(manifest_path).read_text())
    record = next(r for r in manifest["files"] if r["product"] == "mutantcensus")
    source_hash = sha256(archive)
    if source_hash != record["sha256"]:
        raise ValueError("MutantCensus archive SHA256 mismatch; use the pinned v101 GRCh38 TAR.")
    expected_member = record["filename"].replace("_Tsv_", "_").removesuffix(".tar") + ".tsv.gz"
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".cosmic-brca1-", dir=out_dir.parent))
    try:
        counts = Counter()
        transcripts = Counter()
        mutation_types = Counter()
        with tarfile.open(archive, "r:*") as bundle:
            members = bundle.getmembers()
            for member in members:
                path = PurePosixPath(member.name)
                if path.is_absolute() or ".." in path.parts or member.issym() or member.islnk():
                    raise ValueError("Unexpected unsafe TAR member.")
            tables = [m for m in members if m.isfile() and PurePosixPath(m.name).name == expected_member]
            if len(tables) != 1:
                raise ValueError("Expected exactly one " + expected_member + " in the archive.")
            table = tables[0]
            with bundle.extractfile(table) as compressed, gzip.GzipFile(fileobj=compressed) as unzipped:
                with io.TextIOWrapper(unzipped, encoding="utf-8-sig", newline="") as text:
                    reader = csv.DictReader(text, delimiter="\t")
                    fields = reader.fieldnames or []
                    if not fields or len(set(fields)) != len(fields):
                        raise ValueError("Missing or duplicate TSV column names.")
                    gene_field = "GENE_SYMBOL" if "GENE_SYMBOL" in fields else "GENE_NAME"
                    if gene_field not in fields:
                        raise ValueError("No supported gene-symbol column in TSV.")
                    with (staging / "brca1_source_rows.tsv").open("w", newline="") as target:
                        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                        for row in reader:
                            counts["source_rows"] += 1
                            if None in row or any(value is None for value in row.values()):
                                raise ValueError("Malformed TSV row at data row " + str(counts["source_rows"]))
                            if row[gene_field].strip() != "BRCA1":
                                continue
                            writer.writerow(row)
                            counts["brca1_rows"] += 1
                            transcript = row.get("TRANSCRIPT_ACCESSION", row.get("ACCESSION_NUMBER", ""))
                            transcripts[transcript.strip() or "MISSING"] += 1
                            mutation_types[row.get("MUTATION_DESCRIPTION", "").strip() or "MISSING"] += 1
            if counts["brca1_rows"] == 0:
                raise ValueError("No exact BRCA1 gene rows found; inspect the input/schema.")
            archive_members = [{"name": m.name, "bytes": m.size, "regular_file": m.isfile()} for m in members]
        if sha256(archive) != source_hash:
            raise ValueError("Source archive changed during processing.")
        with (staging / "transcript_counts.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["TRANSCRIPT_ACCESSION", "source_rows"])
            writer.writerows(sorted(transcripts.items()))
        report = {
            "stage": "cosmic_brca1_source_extraction", "preparation_passed": True,
            "release": manifest["release"], "genome": manifest["genome"],
            "historical_release_confirmed": False,
            "archive": str(archive.resolve()), "archive_sha256": source_hash,
            "official_archive_hash_matched": True,
            "source_manifest_sha256": sha256(manifest_path), "script_sha256": sha256(__file__),
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_table_member": table.name, "archive_members": archive_members,
            "source_columns": fields, "counts": dict(counts),
            "transcripts": dict(sorted(transcripts.items())),
            "mutation_descriptions": dict(sorted(mutation_types.items())),
            "breast_filter_applied": False, "snv_filter_applied": False,
            "transcript_coordinate_matching_performed": False,
            "observed_cosmic_variant_overlay_reproduced": False,
            "classification_required_for_breast_selection": "PRIMARY_SITE" not in fields,
            "classification_product_url": "https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/classification",
            "interpretation": "All exact-BRCA1 source rows retained, including indels, duplicate observations and every tissue/transcript. These rows are not yet the S12 breast-cancer variant set.",
            "outputs": [{"filename": name, "sha256": sha256(staging / name)}
                        for name in ["brca1_source_rows.tsv", "transcript_counts.csv"]],
        }
        (staging / "cosmic_brca1_preparation_summary.json").write_text(json.dumps(report, indent=2) + "\n")
        if out_dir.exists():
            raise FileExistsError("Output directory appeared during processing.")
        os.rename(staging, out_dir)
        return report
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    report = prepare(args.archive, args.out_dir)
    print(json.dumps({key: report[key] for key in [
        "stage", "preparation_passed", "release", "genome", "official_archive_hash_matched",
        "counts", "transcripts", "breast_filter_applied", "classification_required_for_breast_selection"]}, indent=2))
    print("Report:", args.out_dir / "cosmic_brca1_preparation_summary.json")


if __name__ == "__main__":
    main()
