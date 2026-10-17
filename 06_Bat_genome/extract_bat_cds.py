"""Extract CDS-plus sequences for the frozen historical Bat1K transcript cohorts."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform
import tempfile
import time

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SPECIES = ('molMol2', 'myoMyo6')
HERE = Path(__file__).resolve().parent


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_cohort(path):
    ids = Path(path).read_text().splitlines()
    require(ids and len(ids) == len(set(ids)), 'Empty or duplicate cohort identifiers')
    require(all(name and name.isascii() and all(c.isalnum() or c in '_.-' for c in name)
                for name in ids), 'Unsafe cohort identifier')
    return set(ids)


def read_bed(path, cohort):
    """BED12 uses zero-based, half-open genomic coordinates on both strands."""
    selected, seen = {}, set()
    with Path(path).open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.rstrip('\r\n').split('\t')
            require(len(fields) >= 12, f'Invalid BED12 row {line_number}')
            chrom, start, end, name, _, strand, cds_start, cds_end = fields[:8]
            require(name not in seen, f'Duplicate BED ID: {name}')
            seen.add(name)
            if name not in cohort:
                continue
            start, end, cds_start, cds_end, count = map(int, (start, end, cds_start, cds_end, fields[9]))
            sizes = [int(value) for value in fields[10].rstrip(',').split(',')]
            offsets = [int(value) for value in fields[11].rstrip(',').split(',')]
            require(strand in ('+', '-') and 0 <= start <= cds_start < cds_end <= end,
                    f'Invalid coding interval: {name}')
            require(count > 0 and len(sizes) == len(offsets) == count, f'Invalid blocks: {name}')
            exons = [(start + offset, start + offset + size) for offset, size in zip(offsets, sizes)]
            require(all(start <= a < b <= end for a, b in exons)
                    and all(exons[i][1] <= exons[i + 1][0] for i in range(count - 1)),
                    f'Overlapping/out-of-bounds blocks: {name}')
            require(exons[0][0] == start and exons[-1][1] == end,
                    f'Outer span does not match exon blocks: {name}')
            # Freeze the original study population rather than adding UTR-bearing transcripts.
            require(start == cds_start and end == cds_end,
                    f'Historical cohort has different transcript/CDS bounds: {name}')
            blocks = [(max(a, cds_start), min(b, cds_end)) for a, b in exons
                      if max(a, cds_start) < min(b, cds_end)]
            selected[name] = dict(transcript_id=name, chrom=chrom, strand=strand,
                                  cds_start=cds_start, cds_end=cds_end,
                                  exon_count=count, blocks=blocks)
    require(set(selected) == cohort,
            f'Cohort IDs absent from BED: {sorted(cohort - set(selected))[:10]}')
    return selected


def extract_sequence(genome, row):
    """Return transcript-oriented CDS plus one adjacent genomic base on each end."""
    start, end = row['cds_start'], row['cds_end']
    require(0 < start < end < len(genome), f'Missing flank/outside genome: {row["transcript_id"]}')
    forward = str(genome[start - 1]) + ''.join(str(genome[a:b]) for a, b in row['blocks']) + str(genome[end])
    sequence = forward if row['strand'] == '+' else str(Seq(forward).reverse_complement())
    return sequence.upper()


def inspect_sequence(sequence):
    cds = sequence[1:-1]
    errors = []
    if set(sequence) - set('ACGT'):
        errors.append('non_ACGT')
    if len(cds) < 6 or len(cds) % 3:
        errors.append('invalid_cds_length')
    if cds[:3] != 'ATG':
        errors.append('missing_ATG')
    if cds[-3:] not in ('TAA', 'TAG', 'TGA'):
        errors.append('missing_terminal_stop')
    protein = str(Seq(cds).translate()) if not errors else ''
    return dict(cds_nt=len(cds), protein_aa_including_stop=len(protein),
                internal_stop_count=protein[:-1].count('*'), errors=errors)


def extract(raw_dir, out_dir, species, sources_file=HERE / 'bat1k_pilot_sources.json',
            cohort_manifest_file=HERE / 'cohorts/manifest.json'):
    started = time.perf_counter()
    raw_dir, out_dir = Path(raw_dir).resolve(), Path(out_dir).resolve()
    sources_file, cohort_manifest_file = Path(sources_file).resolve(), Path(cohort_manifest_file).resolve()
    require(species in SPECIES, 'Unknown bat assembly')
    if out_dir.exists():
        raise FileExistsError(f'Choose a new output directory: {out_dir}')
    cohorts = json.loads(cohort_manifest_file.read_text())['cohorts']
    entry = cohorts[species]
    cohort_path = cohort_manifest_file.parent / entry['filename']
    require(cohort_path.parent.resolve() == cohort_manifest_file.parent, 'Invalid cohort path')
    require(digest(cohort_path) == entry['sha256'], 'Cohort SHA256 mismatch')
    cohort = load_cohort(cohort_path)
    require(len(cohort) == entry['transcripts'], 'Cohort count mismatch')
    paths = [raw_dir / f'HL{species}.Bat1Kannotation.bed', raw_dir / f'HL{species}.fa.gz']
    assets = {item['filename']: item for item in json.loads(sources_file.read_text())['assets']}
    downloaded = {}
    download_manifest = raw_dir / 'download_manifest.json'
    if download_manifest.exists():
        manifest = json.loads(download_manifest.read_text())
        require(manifest['download_passed'], 'Raw download manifest did not pass')
        downloaded = {item['filename']: item for item in manifest['outputs']}
    inputs = []
    for path in paths:
        asset = assets[path.name]
        require(path.stat().st_size == asset['expected_bytes'], f'Input size mismatch: {path.name}')
        actual = digest(path)
        require(not asset.get('sha256') or actual == asset['sha256'], f'Input SHA256 mismatch: {path.name}')
        if downloaded:
            require(path.name in downloaded and downloaded[path.name]['sha256'] == actual,
                    f'Download manifest mismatch: {path.name}')
        inputs.append(dict(filename=path.name, sha256=actual, bytes=path.stat().st_size,
                           reference_sha256_checked=bool(asset.get('sha256'))))
    selected = read_bed(paths[0], cohort)
    by_chrom = defaultdict(list)
    for row in selected.values():
        by_chrom[row['chrom']].append(row)
    sequences, rows, seen_contigs = {}, {}, set()
    # Read one contig at a time; no uncompressed whole-genome copy is written.
    with gzip.open(paths[1], 'rt') as stream:
        for genome in SeqIO.parse(stream, 'fasta'):
            require(genome.id not in seen_contigs, f'Duplicate genome contig: {genome.id}')
            seen_contigs.add(genome.id)
            for row in by_chrom.get(genome.id, []):
                name = row['transcript_id']
                try:
                    sequence = extract_sequence(genome.seq, row)
                    qc = inspect_sequence(sequence)
                    if not qc['errors']:
                        sequences[name] = sequence
                except ValueError as error:
                    qc = dict(cds_nt=0, protein_aa_including_stop=0, internal_stop_count=0, errors=[str(error)])
                rows[name] = dict(**{k: v for k, v in row.items() if k != 'blocks'}, **qc)
            if genome.id in by_chrom:
                print(f'{species}: checked {len(rows)}/{len(cohort)} transcripts', flush=True)
    for name in sorted(cohort - set(rows)):
        row = selected[name]
        rows[name] = dict(**{k: v for k, v in row.items() if k != 'blocks'}, cds_nt=0,
                         protein_aa_including_stop=0, internal_stop_count=0, errors=['missing_contig'])
    failed = [row for row in rows.values() if row['errors']]
    report = dict(extraction_passed=not failed, species=species, expected_transcripts=len(cohort),
                  valid_transcripts=len(sequences), failed_transcripts=len(failed),
                  internal_stop_transcripts=sum(row['internal_stop_count'] > 0 for row in rows.values()),
                  internal_stop_policy='retain and flag; historical preprocessing did not exclude them',
                  errors=dict(Counter(error for row in rows.values() for error in row['errors'])),
                  cohort_policy='frozen historical result IDs; no new isoforms or orthologue filtering',
                  flank_policy='one adjacent genomic nucleotide on each side; all bases uppercased',
                  coordinate_policy='BED12 coding-block intersections; negative-strand reverse complement',
                  inputs=inputs, raw_directory=str(raw_dir),
                  sources_sha256=digest(sources_file), cohort_sha256=digest(cohort_path),
                  cohort_manifest_sha256=digest(cohort_manifest_file), script_sha256=digest(__file__),
                  download_manifest_sha256=digest(download_manifest) if download_manifest.exists() else None,
                  python=platform.python_version(), biopython=Bio.__version__,
                  elapsed_seconds=time.perf_counter() - started)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.bat-cds-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'processed'
        stage.mkdir()
        table = stage / 'transcript_qc.csv'
        with table.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(next(iter(rows.values()))))
            writer.writeheader()
            for name in sorted(rows):
                writer.writerow(dict(rows[name], errors=';'.join(rows[name]['errors'])))
        report['qc_sha256'] = digest(table)
        if not failed:
            fasta = stage / f'{species}.cdsplus.fa'
            SeqIO.write((SeqRecord(Seq(sequences[name]), id=name, description='')
                         for name in sorted(sequences)), fasta, 'fasta')
            report.update(fasta_filename=fasta.name, fasta_sha256=digest(fasta))
        (stage / 'extraction_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
        if out_dir.exists():
            raise FileExistsError(f'Output appeared during extraction: {out_dir}')
        stage.rename(out_dir)
    require(not failed, f'{species}: {len(failed)} invalid cohort transcripts; inspect {out_dir}/transcript_qc.csv')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path, help='New species-specific directory')
    parser.add_argument('--species', choices=SPECIES, required=True)
    args = parser.parse_args(argv)
    report = extract(args.raw_dir, args.out_dir, args.species)
    print(json.dumps({key: report[key] for key in ('extraction_passed', 'species', 'valid_transcripts',
                                                    'internal_stop_transcripts', 'elapsed_seconds')}, indent=2))


if __name__ == '__main__':
    main()
