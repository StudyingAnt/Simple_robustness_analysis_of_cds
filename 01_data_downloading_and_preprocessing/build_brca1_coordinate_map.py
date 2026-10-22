#!/usr/bin/env python3
"""Rebuild the small BRCA1 coding-coordinate map from the exact public GENCODE v40 GTF."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re

from download_cosmic_variants import sha256

TRANSCRIPT = "ENST00000357654.9"
CDS_SHA256 = "2a8bfbf3a33c29bc232c0624f2e685fbe1ca2f114d7e452e1e09d87c597e1929"
CDSPLUS_SHA256 = "6f6ab596a2d5b4a38126653a5c0fe7b7dfccd61035e98a3fa74e1c132ea87c16"
GTF_SHA256 = "6254e1c52470d74c37f5a7969ee5b1b7debdd9d6bbd9bf722e65fc3d873f3104"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"
DEFAULT_MAP = Path(__file__).with_name("brca1_gencode_v40_coordinates.json")


def build(gtf):
    if sha256(gtf) != GTF_SHA256:
        raise ValueError("Expected the pinned comprehensive chromosome GENCODE v40 GTF.")
    features = []
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#") or ('transcript_id "' + TRANSCRIPT + '"') not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in ("CDS", "stop_codon"):
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            if attributes.get("gene_name") != "BRCA1" or fields[0] != "chr17" or fields[6] != "-":
                raise ValueError("Unexpected BRCA1 target locus/strand.")
            features.append({"feature": fields[2], "genome_start": int(fields[3]),
                             "genome_end": int(fields[4]), "phase": fields[7]})
    positions = set()
    stop_positions = set()
    for feature in features:
        group = set(range(feature["genome_start"], feature["genome_end"] + 1))
        positions.update(group)
        if feature["feature"] == "stop_codon":
            stop_positions.update(group)
    if len(positions) != 5592 or len(stop_positions) != 3:
        raise ValueError("Target GTF coding length/stop codon does not match the simulated CDS.")
    intervals = []
    for position in sorted(positions):
        if intervals and position == intervals[-1][1] + 1:
            intervals[-1][1] = position
        else:
            intervals.append([position, position])
    offset = 1
    segments = []
    for start, end in reversed(intervals):
        length = end - start + 1
        segments.append({"genome_start": start, "genome_end": end,
                         "cds_start": offset, "cds_end": offset + length - 1})
        offset += length
    return {"schema_version": 1, "assembly": "GRCh38.p13", "gencode_release": 40,
            "gene": "BRCA1", "transcript": TRANSCRIPT, "chromosome": "17", "strand": "-",
            "cds_nt": 5592, "includes_terminal_stop_codon": True,
            "cds_sha256": CDS_SHA256, "cdsplus_sha256": CDSPLUS_SHA256,
            "coordinate_convention": "1-based inclusive genomic and CDS positions; CDS ordered in transcript orientation",
            "source_gtf": {"url": GTF_URL, "sha256": GTF_SHA256,
                           "official_md5": "14a867b82917c8c3006838c3a5053a3e",
                           "md5_manifest_url": GTF_URL.rsplit("/", 1)[0] + "/MD5SUMS"},
            "segments": segments, "source_features": sorted(features, key=lambda f: (f["genome_start"], f["feature"]))}


def validate_genome(mapping, genomic_region):
    """Validate the complete map using an independent fixed-accession genomic FASTA."""
    lines = Path(genomic_region).read_text().splitlines()
    match = re.match(r">NC_000017\.11:(\d+)-(\d+)(?:\s|$)", lines[0])
    if not match or any(line.startswith(">") for line in lines[1:]):
        raise ValueError("Expected one forward-strand NC_000017.11 region FASTA.")
    start, end = int(match[1]), int(match[2])
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if len(sequence) != end - start + 1:
        raise ValueError("Genomic region length mismatch.")
    fragments = []
    for segment in mapping["segments"]:
        if segment["genome_start"] < start or segment["genome_end"] > end:
            raise ValueError("Genomic region does not cover the complete target CDS.")
        fragment = sequence[segment["genome_start"] - start:segment["genome_end"] - start + 1]
        if mapping["strand"] == "-":
            fragment = fragment[::-1].translate(str.maketrans("ACGT", "TGCA"))
        fragments.append(fragment)
    assembled = "".join(fragments)
    digest = hashlib.sha256(assembled.encode()).hexdigest()
    if digest != mapping["cds_sha256"] or len(assembled) != mapping["cds_nt"]:
        raise ValueError("Genomic mapping does not reproduce the pinned target CDS.")
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000017.11"
           "&seq_start=" + str(start) + "&seq_stop=" + str(end) + "&rettype=fasta&retmode=text")
    return {"url": url, "chromosome_accession": "NC_000017.11", "genome_start": start,
            "genome_end": end, "region_nt": len(sequence),
            "region_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "mapped_cds_nt": len(assembled), "mapped_cds_sha256": digest, "pinned_target_cds_match": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--out-file", type=Path, required=True)
    parser.add_argument("--genomic-region", type=Path, help="Optional independent NC_000017.11 region FASTA")
    parser.add_argument("--validation-out", type=Path, help="Save independent validation JSON (requires --genomic-region)")
    args = parser.parse_args()
    result = build(args.gtf)
    if args.validation_out and not args.genomic_region:
        parser.error("--validation-out requires --genomic-region")
    validation = validate_genome(result, args.genomic_region) if args.genomic_region else None
    with args.out_file.open("x") as target:
        target.write(json.dumps(result, indent=2) + "\n")
    if args.validation_out:
        with args.validation_out.open("x") as target:
            target.write(json.dumps(validation, indent=2) + "\n")
    print("Coordinate map:", args.out_file, "coding nt:", result["cds_nt"])


if __name__ == "__main__":
    main()
