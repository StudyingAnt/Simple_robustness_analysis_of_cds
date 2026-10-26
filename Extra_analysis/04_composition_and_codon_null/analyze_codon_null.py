"""Audit and summarize a stratified, protein-preserving codon null experiment."""
import argparse
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'03_baselines_and_stability'))
from common import np, require, digest, write_csv, prior


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--simulation',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--timing-dir',type=Path)
    parser.add_argument('--bootstrap',type=int,default=10000)
    args=parser.parse_args()
    require(not args.out_dir.exists() and args.bootstrap>=2000,'Output exists or too few bootstrap replicates')
    started=time.perf_counter()
    source=json.loads((args.inputs/'input_manifest.json').read_text())
    for entry in source['outputs']:
        require(digest(args.inputs/entry['filename'])==entry['sha256'],'Input hash mismatch')
    report=json.loads((args.simulation/'simulation_summary.json').read_text())
    config=json.loads((args.simulation/'configuration.json').read_text())
    require(report['simulation_passed'] and not report['timing_only'] and not config['timing_only'],'Incomplete/timing run')
    require(digest(args.simulation/'configuration.json')==report['configuration_sha256'],'Configuration changed')
    require(config['input_manifest_sha256']==digest(args.inputs/'input_manifest.json'),'Input provenance mismatch')
    with np.load(args.inputs/'inputs.npz',allow_pickle=False) as d:
        sample=d['null_selected']; ids=d['ids'][sample]; genes=d['genes'][sample]
        strata=d['strata'][sample]; all_strata=d['strata']; profiles=d['profiles'].tolist()
        historical=d['means'][sample]
    require(config['selected_transcripts']==ids.tolist() and config['profiles']==profiles,'Selection/profile mismatch')
    n, r, k, t=len(ids),config['recodings'],len(profiles),config['n_trajectories']
    require(r>=2 and t>=2,'At least two recodings and paths are required for variance estimates')
    byid={str(tid):i for i,tid in enumerate(ids)}
    means=np.full((n,r+1,k),np.nan); variances=np.full_like(means,np.nan); probabilities=np.full_like(means,np.nan)
    one_success=np.zeros(k); one_expect=np.zeros(k); one_variance=np.zeros(k)
    seen=set(); protein_hash={}; replay_checks=0
    for item in report['outputs']:
        path=args.simulation/item['filename']
        require(digest(path)==item['sha256'],'Checkpoint hash mismatch')
        with np.load(path,allow_pickle=False) as saved:
            counts=saved['counts']; meta=json.loads(str(saved['metadata'].item()))
        i,v=byid[meta['Transcript']],meta['variant']
        require((i,v) not in seen and 0<=v<=r,'Unexpected/duplicate variant')
        seen.add((i,v))
        require(counts.shape==(t,k) and np.issubdtype(counts.dtype,np.integer) and (counts>=0).all(),'Invalid counts')
        require(meta['protein_preserved'] and meta['trajectories']==t,'Invalid metadata')
        if i in protein_hash:require(meta['protein_sha256']==protein_hash[i],'Protein hash varies within gene')
        protein_hash[i]=meta['protein_sha256']
        means[i,v]=counts.mean(axis=0);variances[i,v]=counts.var(axis=0,ddof=1)
        p=np.asarray(meta['weighted_one_step']);probabilities[i,v]=p
        require(np.isfinite(p).all() and (p>=0).all() and (p<=1).all(),'Invalid exact probability')
        one_success+=(counts>=1).sum(axis=0);one_expect+=t*p;one_variance+=t*p*(1-p)
        if args.timing_dir and (args.timing_dir/path.name).exists():
            with np.load(args.timing_dir/path.name,allow_pickle=False) as old:
                require(np.array_equal(counts[:len(old['counts'])],old['counts']),'Seed replay failed')
            replay_checks+=1
    require(len(seen)==n*(r+1) and np.isfinite(means).all(),'Missing simulation variants')
    z=(one_success-one_expect)/np.sqrt(one_variance)
    require((np.abs(z)<6).all(),'First-step exact-probability audit failed')
    # Equal numbers selected within strata; weight strata by their full representative-cohort sizes.
    strata_ids=sorted(set(strata.tolist()))
    groups=[np.flatnonzero(strata==s) for s in strata_ids]
    sw=np.array([(all_strata==s).mean() for s in strata_ids])
    weights=np.zeros(n)
    for indices,w in zip(groups,sw):weights[indices]=w/len(indices)
    require(np.isclose(weights.sum(),1),'Invalid design weights')
    original=means[:,0]; null=means[:,1:].mean(axis=1); excess=original-null
    # Bootstrap genes WITHIN length x GC3 strata. Recodings/paths are NOT biological replicates.
    rng=np.random.default_rng(20261026)
    bootstrap=np.zeros((args.bootstrap,k))
    for indices,w in zip(groups,sw):
        draws=rng.choice(indices,size=(args.bootstrap,len(indices)),replace=True)
        bootstrap+=w*excess[draws].mean(axis=1)
    ci=np.quantile(bootstrap,[.025,.975],axis=0)
    # Conditional finite-path plus finite-recoding precision, distinct from gene-resampling intervals.
    conditional_variance=variances[:,0]/t+means[:,1:].var(axis=1,ddof=1)/r
    conditional_se=np.sqrt((weights[:,None]**2*conditional_variance).sum(axis=0))
    summaries=[]
    for j,p in enumerate(profiles):
        summaries.append(dict(Signature=p,n_genes=n,recodings_per_gene=r,trajectories_per_sequence=t,
            stratified_weighted_original_mean=float(weights@original[:,j]),
            stratified_weighted_null_mean=float(weights@null[:,j]),
            mean_excess=float(weights@excess[:,j]),gene_bootstrap_ci_low=float(ci[0,j]),gene_bootstrap_ci_high=float(ci[1,j]),
            conditional_path_and_recoding_mc_se=float(conditional_se[j]),
            unweighted_original_mean=float(original[:,j].mean()),unweighted_null_mean=float(null[:,j].mean()),
            fraction_genes_positive_excess=float((excess[:,j]>0).mean()),
            original_vs_saved_mean_difference=float(weights@(original[:,j]-historical[:,j]))))
    comparisons=[]
    for j,p in enumerate(profiles[1:],1):
        distribution=bootstrap[:,0]-bootstrap[:,j]
        limits=np.quantile(distribution,[.025,.975,.05/(2*(k-1)),1-.05/(2*(k-1))])
        comparisons.append(dict(profile_1=profiles[0],profile_2=p,n_genes=n,
            difference_in_mean_excess=float(weights@(excess[:,0]-excess[:,j])),
            pointwise_ci_low=float(limits[0]),pointwise_ci_high=float(limits[1]),
            bonferroni_ci_low=float(limits[2]),bonferroni_ci_high=float(limits[3])))
    gene_rows=[]
    for i in range(n):
        for j,p in enumerate(profiles):
            gene_rows.append(dict(Transcript=str(ids[i]),gene_id=str(genes[i]),stratum=int(strata[i]),
                design_weight=float(weights[i]),Signature=p,original_mean=float(original[i,j]),
                null_mean=float(null[i,j]),excess=float(excess[i,j]),
                null_between_recoding_sd=float(means[i,1:,j].std(ddof=1)),
                original_mc_se=float(np.sqrt(variances[i,0,j]/t)),saved_n1000_mean=float(historical[i,j]),
                original_one_step_probability=float(probabilities[i,0,j]),
                null_mean_one_step_probability=float(probabilities[i,1:,j].mean())) )
    args.out_dir.mkdir(parents=True)
    write_csv(args.out_dir/'codon_null_summary.csv',summaries)
    write_csv(args.out_dir/'sbs2_excess_comparisons.csv',comparisons)
    write_csv(args.out_dir/'gene_level_null_results.csv',gene_rows)
    fig,axes=prior.paper.plt.subplots(1,2,figsize=(10,4))
    a=np.array([s['stratified_weighted_original_mean'] for s in summaries]);b=np.array([s['stratified_weighted_null_mean'] for s in summaries])
    positions=np.arange(k)
    axes[0].bar(positions-.18,a,.36,label='Original CDS');axes[0].bar(positions+.18,b,.36,label='Synonymous-recoding null')
    axes[0].set(xticks=positions,xticklabels=profiles,ylabel='Stratified-weighted mean RUMS');axes[0].legend(frameon=False,fontsize=8)
    e=np.array([s['mean_excess'] for s in summaries])
    axes[1].errorbar(positions,e,yerr=np.vstack([e-ci[0],ci[1]-e]),fmt='o',capsize=3)
    axes[1].axhline(0,color='gray',linewidth=.8)
    axes[1].set(xticks=positions,xticklabels=profiles,ylabel='Original minus null (95% gene-bootstrap CI)')
    prior.paper.save_plot(fig,args.out_dir,'synonymous_codon_null')
    result=dict(analysis_passed=True,stage='synonymous_codon_recoding_null',genes=n,profiles=profiles,
        recodings_per_gene=r,trajectories_per_sequence=t,total_trajectories=n*(r+1)*k*t,
        complete_paths=n*(r+1)*k*t,incomplete_paths=0,protein_preserving_variants=len(seen),
        first_step_oracle_audit=dict(passed=True,aggregate_z_by_profile=dict(zip(profiles,z.tolist())),threshold_abs_z=6),
        timing_to_production_seed_replay_variants=replay_checks,
        summary=summaries,sbs2_excess_comparisons=comparisons,bootstrap_replicates=args.bootstrap,
        bootstrap_scope='Genes resampled within 16 selection strata; finite simulation estimates treated as the gene-level observations',
        precision_scope='Separate conditional MC SE includes original path sampling and finite null recoding/path sampling; not a biological CI',
        design=source['null_selection']+'; strata weighted by the full representative-cohort proportions',
        interpretation=f'A limited human protein-preserving uniform-codon null; GC, codon usage and context frequencies may change. Differences do not identify selection or establish adaptation; coverage is {n} genes and {k} profiles.',
        historical_results_replaced=False,input_manifest_sha256=digest(args.inputs/'input_manifest.json'),
        simulation_summary_sha256=digest(args.simulation/'simulation_summary.json'),script_sha256=digest(__file__),
        elapsed_seconds=time.perf_counter()-started,
        outputs=[dict(filename=p.name,sha256=digest(p)) for p in sorted(args.out_dir.iterdir())])
    (args.out_dir/'analysis_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='outputs'},indent=2),flush=True)


if __name__=='__main__':main()
