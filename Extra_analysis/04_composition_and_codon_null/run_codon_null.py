"""Protein-preserving synonymous-recoding null using the unchanged Sinabro API."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from Bio.Seq import Seq
from Bio.Data import CodonTable

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'03_baselines_and_stability'))
from common import require, digest, write_csv
from analyze_baselines import signature_table, channel_lookup, opportunities


def seed_for(master, *parts):
    return int.from_bytes(hashlib.sha256(('\0'.join(map(str, (master, *parts)))).encode()).digest()[:16], 'big')


def initialize(source, expected_hash):
    global sinabro, LOOKUP, WEIGHTS, SYNONYMS
    sys.path.insert(0, source)
    import sinabro
    require(sinabro.__version__ == '3.0.0', 'Unexpected Sinabro version')
    package = Path(sinabro.__file__).resolve().parent
    h = hashlib.sha256()
    for path in sorted(package.rglob('*.py')):
        h.update(str(path.relative_to(package.parent)).encode()+b'\0'); h.update(path.read_bytes())
    require(h.hexdigest() == expected_hash, 'Package differs from audited human run')
    labels, WEIGHTS = signature_table(package/'data/COSMIC_v3.3.1_SBS_GRCh37.txt', PROFILES)
    LOOKUP = channel_lookup(labels)
    table = CodonTable.unambiguous_dna_by_id[1]
    SYNONYMS = {}
    for codon in sorted(table.forward_table):
        SYNONYMS.setdefault(table.forward_table[codon], []).append(codon)
    SYNONYMS['*'] = sorted(table.stop_codons)


PROFILES = ['SBS2', 'SBS13', 'SBS19', 'SBS1', 'SBS4']


def compute(job):
    identifier, original, variant, ntraj, master, max_steps = job
    started = time.perf_counter()
    protein = str(Seq(original[1:-1]).translate(table=1, to_stop=False))
    sequence = original
    if variant:
        rng = np.random.default_rng(seed_for(master, 'recoding', identifier, variant))
        sequence = original[0] + ''.join(SYNONYMS[aa][rng.integers(len(SYNONYMS[aa]))] for aa in protein) + original[-1]
    require(len(sequence) == len(original) and sequence[0] == original[0] and sequence[-1] == original[-1], 'Flank/length changed')
    require(str(Seq(sequence[1:-1]).translate(table=1, to_stop=False)) == protein, 'Protein changed on recoding')
    counts = np.empty((ntraj, len(PROFILES)), dtype=np.int32)
    for j, signature in enumerate(PROFILES):
        seeds = np.random.SeedSequence(seed_for(master, 'paths', identifier, variant, signature)).spawn(ntraj)
        for k, seed in enumerate(seeds):
            trajectory = sinabro.Trajectory(identifier, sequence, seed=seed, cds_start=1, cds_end=len(sequence)-1)
            trajectory.autofill(condition='nonsynonymous', method='signature', mutational_signature=signature,
                                cosmic_version='3.3.1', genome_ref='GRCh37', max_steps=max_steps)
            require(trajectory.terminal_reason == 'nonsynonymous' and trajectory.accepted_steps is not None,
                    'Incomplete null trajectory: ' + str(trajectory.terminal_reason))
            counts[k,j] = trajectory.accepted_steps
    total, syn = opportunities(sequence, LOOKUP)
    meta = dict(Transcript=identifier, variant=variant, trajectories=ntraj,
        sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        protein_sha256=hashlib.sha256(protein.encode()).hexdigest(), protein_preserved=True,
        changed_codons=sum(original[p:p+3] != sequence[p:p+3] for p in range(1,len(original)-1,3)),
        cds_nt=len(sequence)-2, seconds=time.perf_counter()-started,
        weighted_one_step=((syn@WEIGHTS)/(total@WEIGHTS)).tolist())
    return counts, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--sinabro-source', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--n-traj', type=int, default=200)
    parser.add_argument('--recodings', type=int, default=20)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--seed', type=int, default=20261025)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--timing-genes', type=int, default=0, help='Length-spaced subset for timing only')
    args = parser.parse_args()
    require(min(args.n_traj, args.recodings, args.workers, args.max_steps)>0, 'Invalid configuration')
    require(not args.out_dir.exists(), 'Output exists; choose a fresh directory')
    source_report = json.loads((args.inputs/'input_manifest.json').read_text())
    for entry in source_report['outputs']:
        require(digest(args.inputs/entry['filename']) == entry['sha256'], 'Prepared input hash mismatch')
    with np.load(args.inputs/'inputs.npz', allow_pickle=False) as data:
        ids, sample = data['ids'], data['null_selected']
        seqs = dict(zip(sample, data['null_sequences']))
        require(data['profiles'].tolist() == PROFILES, 'Profile order mismatch')
        if args.timing_genes:
            order = sorted(sample, key=lambda i: len(seqs[i]))
            sample = np.asarray([order[j] for j in np.linspace(0,len(order)-1,args.timing_genes,dtype=int)])
        jobs = [(str(ids[i]), str(seqs[i]), r, args.n_traj, args.seed, args.max_steps)
                for i in sample for r in range(args.recodings+1)]
    args.out_dir.mkdir(parents=True)
    config = dict(input_manifest_sha256=digest(args.inputs/'input_manifest.json'),
        n_trajectories=args.n_traj, recodings=args.recodings, seed=args.seed, max_steps=args.max_steps,
        profiles=PROFILES, selected_transcripts=[str(ids[i]) for i in sample],
        timing_only=bool(args.timing_genes), null='Uniform independent synonymous codons conditional on the complete translated protein, including stop symbols',
        protein_preserved=True, flanks_preserved=True, codon_usage_gc_and_context_not_preserved=True,
        source_package_sha256=source_report['package_source_sha256'], signature_sha256=source_report['signature_sha256'],
        script_sha256=digest(__file__))
    (args.out_dir/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
    started = time.perf_counter()
    metadata = []
    context = mp.get_context('spawn')
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=context, initializer=initialize,
            initargs=(str(args.sinabro_source.resolve()), source_report['package_source_sha256'])) as executor:
        pending = {executor.submit(compute, job): job for job in jobs}
        for future in as_completed(pending):
            counts, meta = future.result()
            filename = f"{meta['Transcript']}.v{meta['variant']:03d}.npz"
            np.savez_compressed(args.out_dir/filename, counts=counts,
                                metadata=np.array(json.dumps(meta,sort_keys=True)))
            meta['filename'] = filename; meta['sha256'] = digest(args.out_dir/filename)
            metadata.append(meta)
            if len(metadata) % 40 == 0 or len(metadata) == len(jobs):
                print(f'Completed {len(metadata)}/{len(jobs)} sequence variants; elapsed {time.perf_counter()-started:.1f}s', flush=True)
    metadata.sort(key=lambda m:(m['Transcript'], m['variant']))
    report = dict(simulation_passed=True, timing_only=bool(args.timing_genes), genes=len(sample),
        sequence_variants=len(jobs), total_trajectories=len(jobs)*args.n_traj*len(PROFILES),
        incomplete_paths=0, protein_checks_passed=len(metadata), workers=args.workers,
        elapsed_seconds=time.perf_counter()-started, summed_worker_seconds=sum(m['seconds'] for m in metadata),
        configuration_sha256=digest(args.out_dir/'configuration.json'), outputs=metadata)
    (args.out_dir/'simulation_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='outputs'},indent=2),flush=True)


if __name__ == '__main__':
    main()
