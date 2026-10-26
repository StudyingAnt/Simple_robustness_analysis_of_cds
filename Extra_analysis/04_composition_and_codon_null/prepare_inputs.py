"""Authenticate saved human RUMS and prepare sequence-only covariates/selection."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '03_baselines_and_stability'))
from common import np, load_human, digest, require, write_csv
from analyze_baselines import CODONS, ENCODE, signature_table, channel_lookup, opportunities
from Bio import SeqIO

PROFILES = ['SBS2', 'SBS13', 'SBS19', 'SBS1', 'SBS4']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--human-export', type=Path, required=True)
    parser.add_argument('--human-fasta', type=Path, required=True)
    parser.add_argument('--signatures', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--genes-per-stratum', type=int, default=6)
    parser.add_argument('--seed', type=int, default=20261025)
    args = parser.parse_args()
    args.n_traj, args.allow_subset = 1000, False
    require(args.genes_per_stratum > 0 and not args.out_dir.exists(), 'Invalid selection/output')
    started = time.perf_counter()
    data, selected = load_human(args)
    require(digest(args.signatures) == data['manifest']['configuration']['signature_sha256'], 'Signature hash mismatch')
    labels, weights = signature_table(args.signatures, PROFILES)
    lookup = channel_lookup(labels)
    chosen = {data['ids'][i]: j for j, i in enumerate(selected)}
    seqs = [''] * len(selected)
    for record in SeqIO.parse(data['fasta'], 'fasta'):
        tid = record.id.split('|')[0]
        if tid in chosen:
            require(not seqs[chosen[tid]], 'Duplicate selected sequence')
            seqs[chosen[tid]] = str(record.seq).upper()
    require(all(seqs), 'Missing selected sequence')
    names = ['log_cds_nt', 'gc_total', 'gc3'] + ['codon_' + c for c in CODONS] + ['context_' + c for c in CODONS]
    x = np.empty((len(selected), len(names)))
    baseline = np.empty((len(selected), len(PROFILES)))
    rows = []
    for j, (i, sequence) in enumerate(zip(selected, seqs)):
        e = ENCODE[np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)].astype(np.int64)
        require((e < 4).all() and (len(e)-2) % 3 == 0, 'Invalid CDS-plus')
        cds = e[1:-1]
        triplets = cds.reshape(-1, 3)
        codons = triplets[:, 0]*16 + triplets[:, 1]*4 + triplets[:, 2]
        contexts = e[:-2]*16 + e[1:-1]*4 + e[2:]
        x[j] = np.r_[np.log(len(cds)), np.isin(cds, [1,2]).mean(),
                     np.isin(triplets[:,2], [1,2]).mean(),
                     np.bincount(codons, minlength=64)/len(codons),
                     np.bincount(contexts, minlength=64)/len(contexts)]
        total, syn = opportunities(sequence, lookup)
        baseline[j] = (syn @ weights)/(total @ weights)
        a = data['annotations'][i]
        rows.append(dict(Transcript=data['ids'][i], gene_id=a['gene_id'], gene_name=a['gene_name'],
                         cds_nt=len(cds), gc_total=x[j,1], gc3=x[j,2]))
        if (j+1) % 3000 == 0:
            print(f'Sequence covariates {j+1}/{len(selected)}', flush=True)
    require(np.isfinite(x).all() and np.isfinite(baseline).all(), 'Nonfinite covariate')
    # Empirical quartile bins; ties stay in a single bin. Selection never reads outcomes.
    lengths = np.array([len(s)-2 for s in seqs])
    length_edges = np.quantile(lengths, [.25,.5,.75])
    gc_edges = np.quantile(x[:,2], [.25,.5,.75])
    strata = np.searchsorted(length_edges, lengths, side='right')*4 + np.searchsorted(gc_edges, x[:,2], side='right')
    rng = np.random.default_rng(args.seed)
    sample = []
    for s in range(16):
        candidates = np.flatnonzero(strata == s)
        require(len(candidates) >= args.genes_per_stratum, 'Insufficient stratum size')
        sample.extend(sorted(rng.choice(candidates, args.genes_per_stratum, replace=False).tolist()))
    sample = np.asarray(sorted(sample), dtype=np.int64)
    for j, row in enumerate(rows):
        row.update(stratum=int(strata[j]), null_selected=bool(j in set(sample.tolist())))
    args.out_dir.mkdir(parents=True)
    write_csv(args.out_dir/'gene_selection.csv', rows)
    # Deterministic hash folds; every profile uses exactly the same gene partitions.
    folds = np.array([int.from_bytes(hashlib.sha256(r['gene_id'].split('.')[0].encode()).digest()[:8], 'big') % 5 for r in rows])
    np.savez_compressed(args.out_dir/'inputs.npz', x=x, feature_names=np.array(names),
        profiles=np.array(PROFILES), baseline=baseline, folds=folds, strata=strata,
        ids=np.array([r['Transcript'] for r in rows]), genes=np.array([r['gene_id'] for r in rows]),
        null_sequences=np.array([seqs[i] for i in sample]), null_selected=sample,
        means=data['matrix'][selected][:, [data['signatures'].index(p) for p in PROFILES]])
    report = dict(preparation_passed=True, representative_genes=len(selected), null_genes=len(sample),
        null_selection=f'{args.genes_per_stratum} genes per length-quartile x GC3-quartile stratum; no RUMS-based selection',
        selection_seed=args.seed, length_quartiles=length_edges.tolist(), gc3_quartiles=gc_edges.tolist(),
        null_cds_nt_range=[min(rows[j]['cds_nt'] for j in sample), max(rows[j]['cds_nt'] for j in sample)],
        profiles=PROFILES, covariates=names, human_export_manifest_sha256=data['manifest_sha256'],
        fasta_sha256=digest(args.human_fasta), signature_sha256=digest(args.signatures),
        package_source_sha256=data['manifest']['configuration']['package_source_sha256'],
        script_sha256=digest(__file__), elapsed_seconds=time.perf_counter()-started,
        outputs=[dict(filename=p.name, sha256=digest(p)) for p in sorted(args.out_dir.iterdir())])
    (args.out_dir/'input_manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('covariates','outputs')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
