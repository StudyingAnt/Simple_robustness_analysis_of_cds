"""Recompute the reported null estimates/intervals from bundled accepted counts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,default=Path(__file__).resolve().parent/'reported_results')
    args=parser.parse_args()
    folder=args.folder
    with (folder/'SHA256SUMS').open() as f:
        for line in f:
            expected,name=line.strip().split('  ',1)
            if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=expected:
                raise ValueError('Reported file changed: '+name)
    with np.load(folder/'codon_null_path_counts.npz',allow_pickle=False) as d:
        counts=d['counts'];profiles=d['profiles'].tolist();strata=d['strata'];weights=d['design_weights']
        probabilities=d['weighted_one_step']
    assert counts.ndim==4 and (counts>=0).all() and np.isclose(weights.sum(),1)
    n,variants,t,k=counts.shape
    means=counts.mean(axis=2)
    original=means[:,0];null=means[:,1:].mean(axis=1);excess=original-null
    rng=np.random.default_rng(20261026)
    bootstrap=np.zeros((10000,k))
    for s in sorted(set(strata.tolist())):
        indices=np.flatnonzero(strata==s)
        draws=rng.choice(indices,size=(10000,len(indices)),replace=True)
        bootstrap+=weights[indices].sum()*excess[draws].mean(axis=1)
    ci=np.quantile(bootstrap,[.025,.975],axis=0)
    with (folder/'codon_null_summary.csv').open() as f:summary=list(csv.DictReader(f))
    assert [r['Signature'] for r in summary]==profiles
    for j,row in enumerate(summary):
        expected=[weights@original[:,j],weights@null[:,j],weights@excess[:,j],ci[0,j],ci[1,j]]
        saved=[float(row[key]) for key in ('stratified_weighted_original_mean','stratified_weighted_null_mean','mean_excess','gene_bootstrap_ci_low','gene_bootstrap_ci_high')]
        np.testing.assert_allclose(expected,saved,atol=1e-12,rtol=1e-12)
    with (folder/'sbs2_excess_comparisons.csv').open() as f:comparisons=list(csv.DictReader(f))
    for j,row in enumerate(comparisons,1):
        values=bootstrap[:,0]-bootstrap[:,j]
        limits=np.quantile(values,[.025,.975,.05/(2*(k-1)),1-.05/(2*(k-1))])
        np.testing.assert_allclose(limits,[float(row[key]) for key in ('pointwise_ci_low','pointwise_ci_high','bonferroni_ci_low','bonferroni_ci_high')],atol=1e-12,rtol=1e-12)
    successes=(counts>=1).sum(axis=(0,1,2))
    expected=t*probabilities.sum(axis=(0,1))
    variance=t*(probabilities*(1-probabilities)).sum(axis=(0,1))
    z=(successes-expected)/np.sqrt(variance)
    assert (np.abs(z)<6).all()
    print(json.dumps(dict(verification_passed=True,genes=n,profiles=profiles,
        total_trajectories=int(counts.size),recomputed_effects_and_confidence_intervals=True,
        first_step_oracle_z=z.tolist()),indent=2))


if __name__=='__main__':main()
