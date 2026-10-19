"""Shared provenance and input checks for the bat COSMIC production workflow."""
import csv
import json
from pathlib import Path

from Bio import SeqIO

from extract_bat_cds import digest, load_cohort, require

HERE = Path(__file__).resolve().parent
SPECIES = ('molMol2', 'myoMyo6')
SIGNATURES = ['SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS6', 'SBS9', 'SBS13', 'SBS18']


def array_case(index, shards):
    require(shards > 0 and 0 <= index < len(SPECIES) * shards, 'Array index outside both species/shards')
    return SPECIES[index % len(SPECIES)], index // len(SPECIES)


def read_inputs(processed_dir, species, input_path=None, allow_subset=False,
                cohort_manifest=HERE / 'cohorts/manifest.json'):
    """Bind the extraction, frozen population, QC and actual FASTA content."""
    processed_dir = Path(processed_dir).resolve()
    manifest_path = processed_dir / 'extraction_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    require(species in SPECIES and manifest['species'] == species and manifest['extraction_passed'] is True,
            'Wrong species or failed extraction manifest')
    cohort_manifest = Path(cohort_manifest).resolve()
    cohort = json.loads(cohort_manifest.read_text())['cohorts'][species]
    ids_path = cohort_manifest.parent / cohort['filename']
    require(ids_path.parent.resolve() == cohort_manifest.parent, 'Invalid cohort file path')
    require(digest(ids_path) == cohort['sha256'] == manifest['cohort_sha256'], 'Cohort hash mismatch')
    require(digest(cohort_manifest) == manifest['cohort_manifest_sha256'], 'Cohort manifest hash mismatch')
    ids = load_cohort(ids_path)
    require(len(ids) == cohort['transcripts'] == manifest['expected_transcripts'] == manifest['valid_transcripts']
            and manifest['failed_transcripts'] == 0, 'Frozen transcript count mismatch')
    fasta = processed_dir / manifest['fasta_filename']
    require(fasta.parent == processed_dir and digest(fasta) == manifest['fasta_sha256'], 'CDS FASTA hash mismatch')
    qc = processed_dir / 'transcript_qc.csv'
    require(digest(qc) == manifest['qc_sha256'], 'Extraction QC hash mismatch')
    with qc.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == len(ids) and {r['transcript_id'] for r in rows} == ids
            and not any(r['errors'] for r in rows), 'QC transcript set or sequence checks failed')
    annotations = {row['transcript_id']: dict(species=species, **row) for row in rows}
    records = list(SeqIO.parse(fasta, 'fasta'))
    require(len(records) == len(ids) and {r.id for r in records} == ids, 'FASTA transcript set mismatch')
    full = {r.id: r for r in records}
    require(all(r.id == r.description and len(r.seq) - 2 == int(annotations[r.id]['cds_nt'])
                for r in records), 'FASTA identifiers or lengths differ from extraction QC')
    actual = Path(input_path).resolve() if input_path is not None else fasta
    scope = 'full_historical_cohort'
    if actual != fasta:
        require(allow_subset, 'Different input requires explicit --allow-subset for a bounded demonstration')
        records = list(SeqIO.parse(actual, 'fasta'))
        require(records and len(records) == len({r.id for r in records}), 'Empty/duplicate subset input')
        require(all(r.id in full and r.description == full[r.id].description and r.seq == full[r.id].seq
                    for r in records), 'Subset sequence differs from validated full input')
        scope = 'explicit_subset_of_historical_cohort'
    selected = sorted((annotations[r.id] for r in records), key=lambda r: r['transcript_id'])
    return dict(manifest=manifest, manifest_path=manifest_path, input_path=actual,
                input_sha256=digest(actual), annotations=selected, scope=scope,
                full_transcripts=len(ids))


def check_config(config, inputs, n_traj=None):
    require(config['runner'] == 'rums-v3' and config['input_sha256'] == inputs['input_sha256'],
            'Simulation does not belong to this validated bat input')
    require(config['signatures'] == SIGNATURES and config['criterion'] == 'protein_identity'
            and config['cosmic_version'] == '3.3.1' and config['genome_ref'] == 'GRCh37'
            and config['cds_start'] == 1 and config['cds_end'] == 'len(sequence)-1'
            and config['translation_table'] == 1 and config['internal_stops_retained_as_preprocessed'] is True,
            'Unexpected bat COSMIC scientific configuration')
    require(config['n_trajectories'] > 0 and config['shards'] > 0 and config['max_steps'] > 0,
            'Invalid simulation dimensions')
    if n_traj is not None:
        require(config['n_trajectories'] == n_traj, 'Unexpected number of trajectories per combination')


def check_pilot(pilot_dir, processed_root):
    path = Path(pilot_dir).resolve() / 'pilot_report.json'
    pilot = json.loads(path.read_text())
    require(pilot['validation_passed'] is True and pilot['n_trajectories'] >= 100,
            'A passing bat pilot with at least 100 paths per combination is required')
    require(pilot['signatures'] == ['SBS2', 'SBS4', 'SBS13'], 'Unexpected pilot signatures')
    for species in SPECIES:
        row = pilot['species'][species]
        require(row['results_equal'] and row['all_paths_complete'] and row['observed_parallel_execution']
                and row['audits'] == {'serial': True, 'parallel': True}, 'Pilot checks did not pass')
        require(row['extraction_manifest_sha256'] == digest(Path(processed_root) / species / 'extraction_manifest.json'),
                'Pilot used a different extraction manifest')
        require(row['selected_fasta_sha256'] == digest(path.parent / species / 'selected_transcripts.fa'),
                'Pilot selected FASTA changed')
    return path
