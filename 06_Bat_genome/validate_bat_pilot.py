"""Check a small, length-stratified bat COSMIC run with one and multiple workers."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
RUNNER_DIR = HERE.parent / '02_run_sinabro'
sys.path.insert(0, str(RUNNER_DIR))
from validate_parallel import summarize_execution
from validate_results import validate_run

SIGNATURES = ['SBS2', 'SBS4', 'SBS13']


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def select_records(processed_dir, count):
    """One deterministic representative at each equally spaced length quantile."""
    processed_dir = Path(processed_dir).resolve()
    manifest = json.loads((processed_dir / 'extraction_manifest.json').read_text())
    require(manifest['extraction_passed'], 'CDS extraction did not pass')
    fasta = processed_dir / manifest['fasta_filename']
    require(fasta.parent == processed_dir, 'Invalid FASTA path in extraction manifest')
    require(digest(fasta) == manifest['fasta_sha256'], 'Processed FASTA SHA256 mismatch')
    records = list(SeqIO.parse(fasta, 'fasta'))
    require(len(records) == manifest['expected_transcripts'] == manifest['valid_transcripts']
            and len({record.id for record in records}) == len(records), 'Processed transcript set mismatch')
    # Internal stops stay in full input; omit them only from the bounded execution demonstration.
    eligible = [r for r in records if '*' not in str(r.seq[1:-1].translate())[:-1]]
    eligible.sort(key=lambda r: (len(r.seq) - 2, r.id))
    require(1 <= count <= len(eligible), 'Insufficient eligible transcripts for pilot')
    chosen = [eligible[((2 * index + 1) * len(eligible)) // (2 * count)] for index in range(count)]
    return manifest, chosen, dict(total=len(records), eligible=len(eligible),
                                  internal_stop_excluded=len(records) - len(eligible))


def compare(serial, parallel, count, n_traj):
    left = {p.name: p for p in (serial / 'transcripts').glob('*.npz')}
    right = {p.name: p for p in (parallel / 'transcripts').glob('*.npz')}
    require(len(left) == count and set(left) == set(right), 'Pilot checkpoint sets differ')
    checked = 0
    for name in sorted(left):
        with np.load(left[name], allow_pickle=False) as a, np.load(right[name], allow_pickle=False) as b:
            for key in ('accepted_steps', 'terminal_codes'):
                require(a[key].shape == b[key].shape == (n_traj, len(SIGNATURES)), 'Pilot shape mismatch')
                require(np.array_equal(a[key], b[key]), f'Per-path serial/parallel mismatch: {name} {key}')
            ma, mb = [json.loads(str(item['metadata'].item())) for item in (a, b)]
            require(ma['config_sha256'] == mb['config_sha256'] and ma['signatures'] == mb['signatures'] == SIGNATURES,
                    'Pilot scientific configuration mismatch')
            checked += a['accepted_steps'].size
    return checked


def run_pilot(processed_root, out_dir, source, workers=4, n_traj=100, count=12,
              seed=42, max_steps=10000, timeout=3600):
    output, processed_root = Path(out_dir).resolve(), Path(processed_root).resolve()
    require(workers >= 2 and min(n_traj, count, max_steps, timeout) > 0 and seed >= 0,
            'Use at least two workers, positive sizes/cap/timeout, and a nonnegative seed')
    output.mkdir(parents=True, exist_ok=False)
    report = dict(validation_passed=False, signatures=SIGNATURES, n_trajectories=n_traj,
                  records_per_species=count, seed=seed, max_steps=max_steps,
                  pilot_script_sha256=digest(__file__), species={})
    try:
        for species in ('molMol2', 'myoMyo6'):
            manifest, records, population = select_records(processed_root / species, count)
            require(manifest['species'] == species, 'Species manifest mismatch')
            target = output / species
            target.mkdir()
            fasta = target / 'selected_transcripts.fa'
            SeqIO.write(records, fasta, 'fasta')
            with (target / 'selected_transcripts.csv').open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['transcript_id', 'cds_nt', 'length_quantile_midpoint'])
                writer.writerows((r.id, len(r.seq) - 2, (i + .5) / count) for i, r in enumerate(records))
            entry = dict(population=population, selected_fasta_sha256=digest(fasta),
                         extraction_manifest_sha256=digest(processed_root / species / 'extraction_manifest.json'),
                         selected_cds_nt=[len(r.seq) - 2 for r in records], execution={}, audits={})
            report['species'][species] = entry
            for label, allocated in [('serial', 1), ('parallel', workers)]:
                run = target / label
                command = [sys.executable, str(RUNNER_DIR / 'run_rums.py'),
                           '--input', str(fasta), '--out-dir', str(run),
                           '--workers', str(allocated), '--n-traj', str(n_traj),
                           '--seed', str(seed), '--max-steps', str(max_steps),
                           '--expected-records', str(count), '--start-method', 'spawn',
                           '--signatures', *SIGNATURES]
                if source is not None:
                    command.extend(['--sinabro-source', str(Path(source).resolve())])
                print(f'{species} {label}: {count} transcripts x {len(SIGNATURES)} signatures x {n_traj} paths', flush=True)
                started = time.perf_counter()
                subprocess.run(command, check=True, timeout=timeout)
                entry['execution'][label] = dict(summarize_execution(run),
                                                invocation_seconds=time.perf_counter() - started)
                audit = validate_run(run, fasta, target / f'{label}_audit', workers=1, expected_records=count)
                require(audit['validation_passed'], f'{species} {label} archive audit failed')
                entry['audits'][label] = True
            entry['checked_trajectories_per_mode'] = compare(target / 'serial', target / 'parallel', count, n_traj)
            entry['results_equal'] = True
            entry['observed_parallel_execution'] = entry['execution']['parallel']['peak_simultaneous_cases'] > 1
            entry['all_paths_complete'] = all(mode['incomplete_paths'] == 0 for mode in entry['execution'].values())
            entry['elapsed_speedup'] = entry['execution']['serial']['invocation_seconds'] / entry['execution']['parallel']['invocation_seconds']
            require(entry['observed_parallel_execution'] and entry['all_paths_complete'],
                    f'{species}: parallel execution or complete paths not confirmed')
        report['validation_passed'] = True
    except Exception as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        (output / 'pilot_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--processed-root', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path)
    parser.add_argument('--sinabro-source', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--n-traj', type=int, default=100)
    parser.add_argument('--records', type=int, default=12)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--timeout', type=int, default=3600, help='Seconds per simulation invocation')
    args = parser.parse_args(argv)
    report = run_pilot(args.processed_root, args.out_dir, args.sinabro_source, args.workers,
                       args.n_traj, args.records, args.seed, args.max_steps, args.timeout)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
