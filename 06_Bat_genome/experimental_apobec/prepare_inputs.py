"""Bind the already-validated human and bat CDS inputs to an experimental campaign."""
import argparse
import json
from pathlib import Path
import sys

from Bio import SeqIO

from prepare_profiles import digest, load_profiles, require

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from bat_inputs import read_inputs

HUMAN_SHA256 = 'c2f30678cb9f2252bd4717d1f0b8fef9e9575b4be843aa64478c4b1e0422a751'
DATASETS = ('human', 'molMol2', 'myoMyo6')
COUNTS = dict(human=73214, molMol2=31185, myoMyo6=23213)


def prepare(human_input, bat_root, profiles_dir, out_dir):
    human_input, bat_root, profiles_dir, out_dir = [Path(p).resolve() for p in
                                                 (human_input, bat_root, profiles_dir, out_dir)]
    require(not out_dir.exists(), 'Choose a new campaign input directory')
    body = load_profiles(profiles_dir)
    require(digest(human_input) == HUMAN_SHA256, 'Human input differs from the checked GENCODE v40 CDS-plus FASTA')
    records = list(SeqIO.parse(human_input, 'fasta'))
    ids = [r.id.split('|')[0] for r in records]
    require(len(records) == len(set(ids)) == COUNTS['human'], 'Human count/identifier mismatch')
    stops = sum('*' in str(r.seq[1:-1].translate())[:-1] for r in records)
    data = dict(human=dict(input_path=str(human_input), input_sha256=HUMAN_SHA256,
                          transcripts=len(records), internal_stop_transcripts=stops))
    for species in DATASETS[1:]:
        inputs = read_inputs(bat_root / species, species)
        require(len(inputs['annotations']) == COUNTS[species], 'Unexpected bat cohort count')
        data[species] = dict(input_path=str(inputs['input_path']), input_sha256=inputs['input_sha256'],
                             transcripts=COUNTS[species], extraction_manifest_path=str(inputs['manifest_path']),
                             extraction_manifest_sha256=digest(inputs['manifest_path']),
                             internal_stop_transcripts=inputs['manifest']['internal_stop_transcripts'])
    manifest = dict(preparation_passed=True, model='experimental_apobec_NNC_v1', datasets=data,
                    profiles_dir=str(profiles_dir), profiles_sha256=digest(profiles_dir / 'profiles.json'),
                    profile_manifest_sha256=digest(profiles_dir / 'preparation_manifest.json'),
                    profile_ids=[row['id'] for row in body['profiles']],
                    source_workbook_sha256=body['source_sha256'], script_sha256=digest(__file__),
                    combinations=sum(COUNTS.values()) * 6, production_paths_at_n1000=sum(COUNTS.values()) * 6000)
    out_dir.mkdir(parents=True, exist_ok=False)
    (out_dir / 'inputs_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def load_inputs(directory, verify_fasta=False):
    path = Path(directory).resolve() / 'inputs_manifest.json'
    body = json.loads(path.read_text())
    require(body['preparation_passed'] is True and body['model'] == 'experimental_apobec_NNC_v1'
            and set(body['datasets']) == set(DATASETS), 'Campaign inputs did not pass')
    profiles_dir = Path(body['profiles_dir'])
    profiles = load_profiles(profiles_dir)
    require(body['profiles_sha256'] == digest(profiles_dir / 'profiles.json')
            and body['profile_manifest_sha256'] == digest(profiles_dir / 'preparation_manifest.json')
            and body['profile_ids'] == [row['id'] for row in profiles['profiles']], 'Campaign profile mismatch')
    for dataset in DATASETS:
        row = body['datasets'][dataset]
        require(row['transcripts'] == COUNTS[dataset], 'Campaign population mismatch')
        if dataset == 'human':
            require(row['input_sha256'] == HUMAN_SHA256, 'Wrong human cohort')
        else:
            require(digest(row['extraction_manifest_path']) == row['extraction_manifest_sha256'],
                    'Bat extraction manifest changed')
        if verify_fasta:
            require(digest(row['input_path']) == row['input_sha256'], f'{dataset} FASTA changed')
    return body


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--human-input', type=Path, required=True)
    parser.add_argument('--bat-root', type=Path, required=True)
    parser.add_argument('--profiles-dir', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(args.human_input, args.bat_root, args.profiles_dir, args.out_dir), indent=2))


if __name__ == '__main__':
    main()
