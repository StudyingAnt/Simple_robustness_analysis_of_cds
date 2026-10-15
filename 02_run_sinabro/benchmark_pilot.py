"""Measure a small serial Sinabro pilot across CDS lengths and SBS signatures."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

BINS = [('short', 0, 1000), ('medium', 1000, 3000),
        ('long', 3000, 10000), ('very_long', 10000, 10**12)]
SIGNATURES = ['SBS1', 'SBS2', 'SBS13']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def git_head(directory):
    return subprocess.check_output(['git', '-C', str(directory), 'rev-parse', 'HEAD'], text=True).strip()


def case_seed(master_seed, transcript_id, signature):
    key = f'{master_seed}\0{transcript_id}\0{signature}'
    return int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], 'big')


def select_transcripts(input_file):
    from Bio import SeqIO
    candidates = {label: [] for label, _, _ in BINS}
    populations = {label: 0 for label, _, _ in BINS}
    exclusions = {label: 0 for label, _, _ in BINS}
    for record in SeqIO.parse(input_file, 'fasta'):
        text = str(record.seq)
        cds = text[1:-1]
        label = next(name for name, low, high in BINS if low < len(cds) <= high)
        populations[label] += 1
        # Timing ordinary complete ORFs first; retain the full preprocessing dataset.
        internal_stop = any(cds[i:i + 3] in ('TAA', 'TAG', 'TGA')
                            for i in range(0, len(cds) - 3, 3))
        if internal_stop:
            exclusions[label] += 1
            continue
        fields = record.id.split('|')
        candidates[label].append(dict(transcript_id=fields[0], gene_id=fields[1],
                                      transcript_name=fields[4], gene=fields[5],
                                      fasta_id=record.id, cds_nt=len(cds), bin=label))
    selected = []
    used_genes = set()
    for label, _, _ in BINS:
        group = sorted(candidates[label], key=lambda row: (row['cds_nt'], row['transcript_id']))
        if len(group) < 3:
            raise ValueError(f'Not enough eligible transcripts in {label}')
        quantiles = [.1, .5, 1.] if label == 'very_long' else [.1, .5, .9]
        for quantile in quantiles:
            target = group[round(quantile * (len(group) - 1))]['cds_nt']
            available = [row for row in group if row['gene_id'] not in used_genes]
            choice = min(available, key=lambda row: (abs(row['cds_nt'] - target), row['transcript_id']))
            selected.append(dict(choice, target_length_quantile=quantile))
            used_genes.add(choice['gene_id'])
    keep = {row['fasta_id'] for row in selected}
    records = [record for record in SeqIO.parse(input_file, 'fasta') if record.id in keep]
    if len(records) != 12:
        raise ValueError('Selection must contain exactly 12 transcripts')
    return selected, records, populations, exclusions


def worker(args):
    setup_start = time.perf_counter()
    sys.path.insert(0, str(args.sinabro_source.resolve()))
    import numpy as np
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq
    import sinabro
    from sinabro.core import full_signature, signature_path
    row = next(record for record in SeqIO.parse(args.input, 'fasta')
               if record.id.split('|')[0] == args.transcript_id)
    seq = str(row.seq)
    original = str(Seq(seq[1:-1]).translate())
    profile_path = signature_path('3.3.1', 'GRCh37')
    # Cache preparation only: no extra trajectories are simulated.
    full_signature(profile_path, args.signature, .5)
    computer = sinabro.RobustnessComputer(args.transcript_id, seq, seed=args.case_seed,
                                         cds_start=1, cds_end=len(seq) - 1)
    setup_seconds = time.perf_counter() - setup_start
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    trajectories = computer.generate_trajectories(
        args.n_traj, 'signature', 'nonsynonymous', mutational_signature=args.signature,
        cosmic_version='3.3.1', genome_ref='GRCh37', max_steps=args.max_steps)
    simulation_seconds = time.perf_counter() - wall_start
    simulation_cpu_seconds = time.process_time() - cpu_start
    validation_start = time.perf_counter()
    paths = []
    for index, trajectory in enumerate(trajectories):
        records = trajectory.get_records()
        path_hash = hashlib.sha256()
        changed_protein_at = []
        for step, record in enumerate(records):
            sequence = str(record.sequence)
            path_hash.update(sequence.encode() + b'\0')
            if step == 0:
                if sequence != seq:
                    raise AssertionError('Initial record differs from input')
                continue
            previous = str(records[step - 1].sequence)
            changed = [i for i, (a, b) in enumerate(zip(previous, sequence)) if a != b]
            if len(sequence) != len(previous) or len(changed) != 1:
                raise AssertionError('Each step must contain exactly one substitution')
            position = changed[0]
            if not 1 <= position < len(seq) - 1:
                raise AssertionError('Mutation outside CDS')
            if sequence[0] != seq[0] or sequence[-1] != seq[-1]:
                raise AssertionError('Flanking bases changed')
            expected = f'c.{position}{previous[position]}>{sequence[position]}'
            if record.hgvs_mrna != expected:
                raise AssertionError('Incorrect CDS-relative mutation coordinate')
            if str(Seq(sequence[1:-1]).translate()) != original:
                changed_protein_at.append(step)
        complete = trajectory.terminal_reason == 'nonsynonymous'
        if complete:
            if changed_protein_at != [trajectory.get_length()]:
                raise AssertionError('Protein identity did not hold before the terminal event')
            if trajectory.accepted_steps != trajectory.get_length() - 1:
                raise AssertionError('Incorrect preterminal count')
            if trajectory.synonymous_steps != trajectory.accepted_steps:
                raise AssertionError('Synonymous and accepted counts differ')
        elif changed_protein_at:
            raise AssertionError('Incomplete path contains a protein change')
        paths.append(dict(trajectory_id=index, terminal_reason=trajectory.terminal_reason,
                          n_steps=trajectory.get_length(), accepted_steps=trajectory.accepted_steps,
                          trajectory_sha256=path_hash.hexdigest()))
    validation_seconds = time.perf_counter() - validation_start
    accepted = [path['accepted_steps'] for path in paths if path['accepted_steps'] is not None]
    n_complete = len(accepted)
    sd = float(np.std(accepted, ddof=1)) if n_complete > 1 else None
    result = dict(transcript_id=args.transcript_id, signature=args.signature,
                  cds_nt=len(seq) - 2, n_traj=args.n_traj, case_seed=args.case_seed,
                  setup_seconds=setup_seconds, simulation_seconds=simulation_seconds,
                  simulation_cpu_seconds=simulation_cpu_seconds,
                  validation_seconds=validation_seconds,
                  seconds_per_trajectory=simulation_seconds / args.n_traj,
                  n_complete=n_complete, n_incomplete=args.n_traj - n_complete,
                  total_mutation_steps=sum(path['n_steps'] for path in paths),
                  max_mutation_steps=max(path['n_steps'] for path in paths),
                  # A partial-path average is deliberately not reported as RUMS.
                  mean_accepted=float(np.mean(accepted)) if n_complete == args.n_traj else None,
                  accepted_sd=sd if n_complete == args.n_traj else None,
                  accepted_mc_se=sd / np.sqrt(n_complete) if n_complete == args.n_traj and sd is not None else None,
                  peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                  python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__,
                  sinabro=sinabro.__version__, sinabro_module=str(Path(sinabro.__file__).resolve()),
                  profile_sha256=sha256(profile_path), paths=paths, all_path_checks_passed=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')


def run(args):
    from Bio import SeqIO
    import numpy as np
    import Bio
    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    selected, records, populations, exclusions = select_transcripts(args.input)
    selected_file = output / 'selected_transcripts.fa'
    SeqIO.write(records, selected_file, 'fasta')
    manifest = dict(input_path=str(args.input.resolve()), input_sha256=sha256(args.input),
                    selected_fasta_sha256=sha256(selected_file), selection=selected,
                    population_by_bin=populations, internal_stop_excluded_from_pilot_by_bin=exclusions,
                    selection_rule='10th/50th/90th length percentile per fixed bin; last bin includes maximum; distinct genes',
                    bin_intervals=[dict(name=name, lower_exclusive=low, upper_inclusive=high)
                                   for name, low, high in BINS],
                    signatures=SIGNATURES, n_traj_per_case=args.n_traj, planned_cases=36,
                    master_seed=args.seed, max_steps=args.max_steps, serial_workers=1,
                    python=platform.python_version(), numpy=np.__version__, biopython=Bio.__version__,
                    platform=platform.platform(), logical_cpus=os.cpu_count(),
                    package_source=str(args.sinabro_source.resolve()),
                    package_commit=git_head(args.sinabro_source),
                    analysis_commit=git_head(Path(__file__).parent),
                    script_sha256=sha256(__file__), command=sys.argv,
                    started_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    trajectory_warmup_count=0, repeats=1)
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    with (output / 'selected_transcripts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    jobs = [(row, signature) for row in selected for signature in SIGNATURES]
    order = np.random.default_rng(args.seed).permutation(len(jobs))
    summaries, failures = [], []
    worker_dir = output / 'cases'
    worker_dir.mkdir()
    env = dict(os.environ)
    for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[key] = '1'
    for position, job_index in enumerate(order, 1):
        row, signature = jobs[int(job_index)]
        identifier = row['transcript_id']
        case_output = worker_dir / f'{identifier}_{signature}.json'
        command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                   '--input', str(selected_file), '--output', str(case_output),
                   '--sinabro-source', str(args.sinabro_source.resolve()),
                   '--transcript-id', identifier, '--signature', signature,
                   '--case-seed', str(case_seed(args.seed, identifier, signature)),
                   '--n-traj', str(args.n_traj), '--max-steps', str(args.max_steps)]
        wall_start = time.perf_counter()
        try:
            completed = subprocess.run(command, env=env, text=True, capture_output=True,
                                       timeout=args.case_timeout)
            wall = time.perf_counter() - wall_start
            if completed.returncode:
                raise RuntimeError(completed.stderr or completed.stdout)
            result = json.loads(case_output.read_text())
            paths = result.pop('paths')
            result.update(gene=row['gene'], transcript_name=row['transcript_name'], bin=row['bin'],
                          subprocess_wall_seconds=wall)
            summaries.append(result)
            print(f"{position:02d}/36 {row['gene']} {row['cds_nt']} nt {signature}: "
                  f"{result['simulation_seconds']:.3f}s simulation; "
                  f"{result['n_complete']}/{args.n_traj} complete", flush=True)
        except (RuntimeError, subprocess.TimeoutExpired) as exc:
            failures.append(dict(transcript_id=identifier, signature=signature,
                                 error=str(exc), elapsed_seconds=time.perf_counter() - wall_start))
            print(f'{position:02d}/36 ERROR {identifier} {signature}: {exc}', flush=True)
        (output / 'summary.json').write_text(json.dumps(summaries, indent=2) + '\n')
        (output / 'failures.json').write_text(json.dumps(failures, indent=2) + '\n')
    if summaries:
        with (output / 'summary.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
            writer.writeheader()
            writer.writerows(summaries)
    manifest.update(finished_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    total_driver_wall_seconds=time.perf_counter() - start,
                    successful_cases=len(summaries), failed_cases=len(failures),
                    trajectories_completed=sum(row['n_complete'] for row in summaries),
                    trajectories_incomplete=sum(row['n_incomplete'] for row in summaries))
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if failures:
        raise SystemExit(f'{len(failures)} cases failed; inspect failures.json')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--sinabro-source', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path)
    parser.add_argument('--n-traj', type=int, default=20)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--case-timeout', type=float, default=180)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--transcript-id', help=argparse.SUPPRESS)
    parser.add_argument('--signature', help=argparse.SUPPRESS)
    parser.add_argument('--case-seed', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--output', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.n_traj < 1 or args.max_steps < 1 or args.seed < 0 or args.case_timeout <= 0:
        parser.error('counts and timeout must be positive; seed must be nonnegative')
    if not args.worker and args.out_dir is None:
        parser.error('--out-dir is required')
    if args.worker:
        worker(args)
    else:
        run(args)


if __name__ == '__main__':
    main()
