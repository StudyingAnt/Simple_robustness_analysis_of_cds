"""Extract CDS-plus with the historical length/start/stop/PAR_Y filters."""
import argparse
from collections import Counter
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import tempfile

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FILTER_BITS = {'incomplete_codon': 1, 'missing_start_codon': 2,
               'missing_terminal_stop': 4, 'par_y_duplicate': 8}


def classify(record):
    """Return the original filter mask and GENCODE 1-based inclusive CDS bounds."""
    annotations = [field for field in record.id.split('|') if field.startswith('CDS:')]
    if len(annotations) != 1:
        raise ValueError(f'Expected one CDS annotation: {record.id}')
    match = re.fullmatch(r'CDS:(\d+)-(\d+)', annotations[0])
    if match is None:
        raise ValueError(f'Invalid CDS annotation: {record.id}')
    start, end = map(int, match.groups())
    if not 1 <= start <= end <= len(record.seq):
        raise ValueError(f'CDS coordinates outside transcript: {record.id}')
    cds = record.seq[start - 1:end]
    mask = 0
    if len(cds) % 3:
        mask |= FILTER_BITS['incomplete_codon']
    if cds[:3] != 'ATG':
        mask |= FILTER_BITS['missing_start_codon']
    if cds[-3:] not in ('TAA', 'TAG', 'TGA'):
        mask |= FILTER_BITS['missing_terminal_stop']
    if record.id.split('|')[0].endswith('PAR_Y'):
        mask |= FILTER_BITS['par_y_duplicate']
    return mask, start, end


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def preprocess(input_file, output_dir, release):
    input_file, output_dir = Path(input_file).resolve(), Path(output_dir).resolve()
    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    if output_dir.exists():
        raise FileExistsError(f'Refusing to overwrite output directory: {output_dir}')
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    prefix = f'gencode.v{release}.pc_transcripts'
    accepted_name = prefix + '.nopary.cdsplus.fa'
    rejected_names = {i: prefix + f'.nopass.parity{i}.fa' for i in range(1, 16)}
    missing_name = prefix + '.nopass.missing_flank.fa'
    counts = Counter()
    reasons = Counter()
    non_acgt = 0
    input_records = 0
    script = Path(__file__).resolve()
    commit = subprocess.check_output(['git', '-C', str(script.parent), 'rev-parse', 'HEAD'], text=True).strip()

    with tempfile.TemporaryDirectory(prefix='.cdsplus-', dir=output_dir.parent) as directory:
        stage = Path(directory)
        with ExitStack() as stack:
            handles = {name: stack.enter_context((stage / name).open('x'))
                       for name in [accepted_name, missing_name, *rejected_names.values()]}
            for record in SeqIO.parse(input_file, 'fasta'):
                input_records += 1
                mask, start, end = classify(record)
                if mask:
                    name = rejected_names[mask]
                    for reason, bit in FILTER_BITS.items():
                        if mask & bit:
                            reasons[reason] += 1
                    SeqIO.write(record, handles[name], 'fasta')
                elif start - 2 < 0 or end + 1 > len(record.seq):
                    # The old script silently omitted these otherwise-valid transcripts.
                    name = missing_name
                    SeqIO.write(record, handles[name], 'fasta')
                else:
                    name = accepted_name
                    sequence = record.seq[start - 2:end + 1]
                    non_acgt += bool(set(str(sequence)) - set('ACGT'))
                    result = SeqRecord(sequence, id=record.id, description='')
                    SeqIO.write(result, handles[name], 'fasta')
                counts[name] += 1

        if input_records == 0:
            raise ValueError('Input FASTA has no records')
        if sum(counts.values()) != input_records:
            raise RuntimeError('Record accounting failed')
        outputs = [{'filename': name, 'fasta_records': counts[name],
                    'bytes': (stage / name).stat().st_size, 'sha256': sha256(stage / name)}
                   for name in [accepted_name, missing_name, *rejected_names.values()]]
        report = {'gencode_release': release, 'input_path': str(input_file),
                  'input_sha256': sha256(input_file), 'input_records': input_records,
                  'accepted_records': counts[accepted_name],
                  'rejected_by_original_filters': sum(counts[name] for name in rejected_names.values()),
                  'otherwise_valid_but_missing_flank': counts[missing_name],
                  'non_acgt_accepted_records': non_acgt, 'filter_bits': FILTER_BITS,
                  'individual_filter_counts_may_overlap': dict(reasons),
                  'outputs': outputs, 'analysis_commit': commit,
                  'script_sha256': sha256(script), 'completed_at_utc': datetime.now(timezone.utc).isoformat(),
                  'python': platform.python_version(), 'biopython': Bio.__version__}
        (stage / 'preprocessing_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
        stage.rename(output_dir)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True, help='Uncompressed GENCODE transcript FASTA')
    parser.add_argument('--out-dir', type=Path, required=True, help='New directory for CDS-plus, rejected records and manifest')
    parser.add_argument('--release', type=int, default=40)
    args = parser.parse_args(argv)
    if args.release < 1:
        parser.error('--release must be positive')
    try:
        report = preprocess(args.input, args.out_dir, args.release)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Input: {report['input_records']}; accepted CDS-plus: {report['accepted_records']}; "
          f"filter rejects: {report['rejected_by_original_filters']}; "
          f"missing flanks: {report['otherwise_valid_but_missing_flank']}")
    print(f'Manifest: {args.out_dir / "preprocessing_manifest.json"}')


if __name__ == '__main__':
    main()
