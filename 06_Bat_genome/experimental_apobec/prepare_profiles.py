"""Convert the original six-enzyme XLSX into per-opportunity NNC mutation weights."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET
import zipfile

WORKBOOK_SHA256 = 'af2e8dad39b7ebb3e36d4e55d388c7550967610398b7361dc06cdddaae5add5a'
ENZYMES = [
    ('B7', 'HLmolMol2_APOBEC3A', 'molMol2', 'D', 'E', 'F'),
    ('B8', 'HLmolMol2_APOBEC3C', 'molMol2', 'G', 'H', 'I'),
    ('B24', 'HLmolMol2_APOBEC1', 'molMol2', 'J', 'K', 'L'),
    ('B22', 'HLmyoMyo6_APOBEC3A3', 'myoMyo6', 'M', 'N', 'O'),
    ('B23', 'HLmyoMyo6_APOBEC3C', 'myoMyo6', 'P', 'Q', 'R'),
    ('B27', 'HLmyoMyo6_APOBEC1', 'myoMyo6', 'S', 'T', 'U'),
]
MOTIFS = [''.join(pair) + 'C' for pair in itertools.product('ACGT', repeat=2)]
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def cells(archive, name, strings):
    root = ET.fromstring(archive.read(name))
    require(not root.findall('.//s:f', NS), 'Workbook formulas are not accepted; preserve original fixed values')
    result = {}
    for cell in root.findall('.//s:c', NS):
        value = cell.find('s:v', NS)
        kind = cell.get('t')
        if kind == 'inlineStr':
            result[cell.get('r')] = ''.join(t.text or '' for t in cell.findall('.//s:t', NS))
        elif value is not None and value.text is not None:
            result[cell.get('r')] = strings[int(value.text)] if kind == 's' else float(value.text)
    return result


def read_workbook(path):
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, 'Corrupt XLSX archive')
        shared = ET.fromstring(archive.read('xl/sharedStrings.xml')) if 'xl/sharedStrings.xml' in archive.namelist() else None
        strings = [''.join(t.text or '' for t in item.findall('.//s:t', NS))
                   for item in shared.findall('s:si', NS)] if shared is not None else []
        data = cells(archive, 'xl/worksheets/sheet1.xml', strings)
        key = cells(archive, 'xl/worksheets/sheet2.xml', strings)
    mapping = {str(key.get(f'A{row}', '')).strip(): str(key.get(f'B{row}', '')).strip() for row in range(3, 9)}
    require(mapping == {sample: gene for sample, gene, *_ in ENZYMES}, 'Unexpected enzyme/sample identities')
    require('average' in str(data.get('C23', '')).lower(), 'Source must identify motif values as averages')
    rows = []
    for number in range(4, 20):
        motif = data.get(f'B{number}')
        require(motif in MOTIFS, f'Invalid NNC motif at source row {number}')
        occurrences = data[f'C{number}']
        require(math.isfinite(occurrences) and occurrences > 0 and int(occurrences) == occurrences,
                'Invalid substrate motif count')
        for sample, gene, species, count_col, total_col, percent_col in ENZYMES:
            count, total, percent = [data[f'{col}{number}'] for col in (count_col, total_col, percent_col)]
            require(all(math.isfinite(v) for v in (count, total, percent))
                    and total > 0 and 0 <= count <= total and 0 <= percent <= 100,
                    'Invalid source counts/percentages')
            require(math.isclose(percent, 100 * count / total, rel_tol=1e-10, abs_tol=1e-12),
                    'Percentage differs from mean count / sequenced molecules')
            rows.append(dict(profile_id=gene, sample_code=sample, enzyme_species=species,
                             motif=motif, substrate_motif_occurrences=int(occurrences),
                             mean_ct_count=count, total_molecules=total, percentage=percent))
    return rows


def prepare(workbook, out_dir, expected_sha=WORKBOOK_SHA256):
    workbook, out_dir = Path(workbook).resolve(), Path(out_dir).resolve()
    require(not out_dir.exists(), 'Choose a new output directory')
    source_sha = digest(workbook)
    require(expected_sha is None or source_sha == expected_sha, 'Original workbook SHA256 mismatch')
    rows = read_workbook(workbook)
    profiles = []
    for sample, gene, species, *_ in ENZYMES:
        group = [row for row in rows if row['profile_id'] == gene]
        require(len(group) == 16 and {r['motif'] for r in group} == set(MOTIFS), 'Missing/duplicate source motif')
        require(len({r['total_molecules'] for r in group}) == 1, 'Within-enzyme molecule totals differ')
        ordered = {row['motif']: row for row in group}
        denominator = math.fsum(row['percentage'] for row in group)
        require(denominator > 0, 'All-zero enzyme activity')
        weights = [ordered[motif]['percentage'] / denominator for motif in MOTIFS]
        profiles.append(dict(id=gene, sample_code=sample, enzyme_species=species, weights=weights))
        for row in group:
            row['relative_per_site_weight'] = row['percentage'] / denominator
    body = dict(format='experimental-apobec-nnc-v1', source_sha256=source_sha,
                motifs=[motif[:2] + '[C>T]' for motif in MOTIFS], profiles=profiles,
                weighting='motif-average percentage / sum of 16 percentages; one weight per eligible opportunity',
                substrate_motif_counts_multiplied=False, strand='both', strand_bias=0.5,
                context_policy='two upstream bases for C>T, reverse complement for G>A; full context must exist',
                flank_policy='use original CDS-plus one-base flanks; never mutate outside the CDS',
                historical_probability_csv_recovered=False)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.apobec-prepare-', dir=out_dir.parent) as temporary:
        stage = Path(temporary) / 'profiles'
        stage.mkdir()
        (stage / 'profiles.json').write_text(json.dumps(body, indent=2) + '\n')
        with (stage / 'source_values.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        report = dict(preparation_passed=True, profiles=6, motifs_per_profile=16,
                      workbook_sha256=source_sha, script_sha256=digest(__file__),
                      historical_probability_csv_recovered=False,
                      outputs=[dict(filename=p.name, sha256=digest(p)) for p in sorted(stage.iterdir())])
        (stage / 'preparation_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
        require(not out_dir.exists(), 'Output appeared during preparation')
        stage.rename(out_dir)
    return report


def load_profiles(directory):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / 'preparation_manifest.json').read_text())
    require(manifest['preparation_passed'] and manifest['workbook_sha256'] == WORKBOOK_SHA256,
            'Profile preparation did not pass for the original workbook')
    hashes = {r['filename']: r['sha256'] for r in manifest['outputs']}
    require(set(hashes) == {'profiles.json', 'source_values.csv'}, 'Unexpected prepared files')
    require(all(digest(directory / name) == sha for name, sha in hashes.items()), 'Prepared profile hash mismatch')
    body = json.loads((directory / 'profiles.json').read_text())
    require(body['format'] == 'experimental-apobec-nnc-v1' and body['source_sha256'] == WORKBOOK_SHA256,
            'Wrong profile format/source')
    require(body['motifs'] == [m[:2] + '[C>T]' for m in MOTIFS]
            and body['substrate_motif_counts_multiplied'] is False
            and body['strand'] == 'both' and body['strand_bias'] == .5, 'Unexpected biological profile policy')
    require([p['id'] for p in body['profiles']] == [row[1] for row in ENZYMES], 'Wrong profile identities/order')
    for profile, enzyme in zip(body['profiles'], ENZYMES):
        weights = profile['weights']
        require(profile['sample_code'] == enzyme[0] and profile['enzyme_species'] == enzyme[2], 'Sample identity mismatch')
        require(len(weights) == 16 and all(math.isfinite(v) and v >= 0 for v in weights)
                and math.isclose(math.fsum(weights), 1., rel_tol=1e-12), 'Invalid normalized weights')
    return body


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(args.workbook, args.out_dir), indent=2))


if __name__ == '__main__':
    main()
