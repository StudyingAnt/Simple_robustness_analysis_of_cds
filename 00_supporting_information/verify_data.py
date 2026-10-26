"""Check the published numerical bundle without rerunning any simulations."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from xml.etree import ElementTree as E
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = json.loads((ROOT / 'data_manifest.json').read_text())
    for record in manifest['files']:
        path = ROOT / record['file']
        data = path.read_bytes()
        require(hashlib.sha256(data).hexdigest() == record['sha256'], str(path))
        plain = gzip.decompress(data) if path.suffix == '.gz' else data
        require(hashlib.sha256(plain).hexdigest() == record['uncompressed_sha256'], 'Decoded hash: ' + str(path))
        require(len(plain) == record['uncompressed_bytes'], 'Decoded size: ' + str(path))
    tables = {}
    for record in manifest['datasets']:
        with gzip.open(ROOT / record['file'], 'rt', newline='') as stream:
            rows = list(csv.reader(stream))
        gene = record['row_unit'] != 'transcript'
        metadata = 2 if gene else 1
        header = rows[0]
        require(len(header) == record['columns'] + metadata, 'Column count')
        require(len(rows) - 1 == record['rows'], 'Row count')
        require(len({r[0] for r in rows[1:]}) == record['rows'], 'Repeated row identifiers')
        require(all(len(r) == len(header) and all(math.isfinite(float(v)) and float(v) >= 0 for v in r[metadata:]) for r in rows[1:]), 'Invalid RUMS value')
        tables[(record['species'], record['process'], gene)] = rows
    for process in ('cosmic', 'apobec'):
        transcript = {r[0]: r[1:] for r in tables['human', process, False][1:]}
        require(all(r[2:] == transcript[r[1]] for r in tables['human', process, True][1:]), 'Gene values differ from representative transcript')
        require(len({r[0].split('.')[0] for r in tables['human', process, True][1:]}) == 19309, 'Stable gene count')
    selections = [{r[0]: r[1] for r in tables['human', p, True][1:]} for p in ('cosmic', 'apobec')]
    require(selections[0] == selections[1], 'Gene representatives differ across processes')
    for species in ('human', 'molMol2', 'myoMyo6'):
        a = {r[0] for r in tables[species, 'cosmic', False][1:]}
        b = {r[0] for r in tables[species, 'apobec', False][1:]}
        require(a == b, 'Transcript cohorts differ across processes')
        with gzip.open(ROOT / 'data_matrices' / (species + '_transcript_annotations.csv.gz'), 'rt', newline='') as stream:
            require(a == {r['transcript_id'] for r in csv.DictReader(stream)}, 'Annotation identifiers differ')
    source = ROOT / 'apobec_preferences'
    profiles = json.loads((source / 'profiles.json').read_text())
    with (source / 'source_values.csv').open(newline='') as stream:
        values = {(r['profile_id'], r['motif']): r for r in csv.DictReader(stream)}
    require(len(values) == 96, 'Expected six profiles with 16 motifs')
    for profile in profiles['profiles']:
        require(len(profile['weights']) == 16 and math.isclose(sum(profile['weights']), 1, abs_tol=1e-12), 'Invalid profile weights')
        for motif, weight in zip(profiles['motifs'], profile['weights']):
            require(math.isclose(float(values[profile['id'], motif[:2] + 'C']['relative_per_site_weight']), weight, abs_tol=1e-15), 'Weight differs from source')
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    def sheet(number, index):
        with ZipFile(ROOT / ('Table_S%d.xlsx' % number)) as archive:
            return archive.read('xl/worksheets/sheet%d.xml' % index)
    require(sheet(3, 4) == sheet(4, 4), 'S3/S4 shared SBS2 rows differ')
    require(len(E.fromstring(sheet(1, 2)).findall('s:sheetData/s:row', ns)) == 3082, 'Expected 3081 human paired comparisons')
    print(json.dumps(dict(verification_passed=True, matrices=len(manifest['datasets']),
        supplementary_tables=6, profiles=6, motifs_per_profile=16,
        shared_sbs2_rows_identical=True, new_trajectories=False), indent=2))


if __name__ == '__main__':
    main()
