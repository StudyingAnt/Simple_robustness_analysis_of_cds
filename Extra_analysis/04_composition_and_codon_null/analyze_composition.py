"""Joint sequence-composition adjustment of saved RUMS, with gene-held-out checks."""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'03_baselines_and_stability'))
from common import np, stats, require, digest, write_csv, prior
from scipy.linalg import lstsq


def fit_predict(x, y, train, test):
    center, scale = x[train].mean(axis=0), x[train].std(axis=0)
    keep = scale > 1e-12
    a = np.column_stack([np.ones(len(train)), (x[train][:,keep]-center[keep])/scale[keep]])
    b = np.column_stack([np.ones(len(test)), (x[test][:,keep]-center[keep])/scale[keep]])
    coefficient, _, rank, _ = lstsq(a, y[train], cond=1e-10, lapack_driver='gelsd')
    return b @ coefficient, int(rank)


def r2(actual, predicted):
    return float(1-np.square(actual-predicted).sum()/np.square(actual-actual.mean()).sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    require(not args.out_dir.exists(), 'Output exists')
    started = time.perf_counter()
    manifest = json.loads((args.inputs/'input_manifest.json').read_text())
    for entry in manifest['outputs']:
        require(digest(args.inputs/entry['filename']) == entry['sha256'], 'Input hash mismatch')
    with np.load(args.inputs/'inputs.npz', allow_pickle=False) as d:
        x, means, p0, folds = d['x'], d['means'], d['baseline'], d['folds']
        profiles, ids, genes, names = d['profiles'].tolist(), d['ids'], d['genes'], d['feature_names'].tolist()
    blocks = [('length_gc_gc3',3,False), ('plus_codon_usage',67,False),
              ('plus_trinucleotide_opportunity',131,False), ('plus_weighted_synonymous_probability',131,True)]
    args.out_dir.mkdir(parents=True)
    summary, prediction_rows = [], []
    log_predictions = {}
    for transform in ('log1p', 'raw'):
        y = np.log1p(means) if transform == 'log1p' else means
        for name, stop, weighted in blocks:
            predicted = np.empty_like(y)
            ranks = [[] for _ in profiles]
            for fold in range(5):
                train, test = np.flatnonzero(folds!=fold), np.flatnonzero(folds==fold)
                if weighted:
                    for j in range(len(profiles)):
                        z = np.column_stack([x[:,:stop], p0[:,j]])
                        predicted[test,j], rank = fit_predict(z, y[:,j], train, test)
                        ranks[j].append(rank)
                else:
                    predicted[test], rank = fit_predict(x[:,:stop], y, train, test)
                    for r in ranks: r.append(rank)
            require(np.isfinite(predicted).all(), 'Nonfinite held-out prediction')
            if transform == 'log1p': log_predictions[name] = predicted
            for j, profile in enumerate(profiles):
                z = np.column_stack([x[:,:stop], p0[:,j]]) if weighted else x[:,:stop]
                fitted, full_rank = fit_predict(z,y[:,j],np.arange(len(y)),np.arange(len(y)))
                summary.append(dict(Signature=profile, transform=transform, model=name, n_genes=len(y),
                    feature_columns=z.shape[1], effective_rank_including_intercept=full_rank,
                    min_training_rank=min(ranks[j]), max_training_rank=max(ranks[j]),
                    training_r_squared=r2(y[:,j],fitted), gene_held_out_r_squared=r2(y[:,j],predicted[:,j]),
                    gene_held_out_rmse=float(np.sqrt(np.square(y[:,j]-predicted[:,j]).mean()))))
                if transform == 'log1p':
                    prediction_rows.extend(dict(Transcript=str(tid),gene_id=str(gene),Signature=profile,
                        model=name, fold=int(fold), observed_log1p_rums=float(actual), predicted_log1p_rums=float(pred))
                        for tid,gene,fold,actual,pred in zip(ids,genes,folds,y[:,j],predicted[:,j]))
            print(f'Completed {transform}: {name}',flush=True)
    # Explicit partial-rank association: residualize BOTH ranked variables on ranked covariates.
    ranked_x = np.column_stack([stats.rankdata(x[:,k]) for k in range(x.shape[1])])
    partial = []
    for j, profile in enumerate(profiles):
        ranked_pair = np.column_stack([stats.rankdata(means[:,j]),stats.rankdata(p0[:,j])])
        fitted, rank = fit_predict(ranked_x,ranked_pair,np.arange(len(x)),np.arange(len(x)))
        residual = ranked_pair-fitted
        partial.append(dict(Signature=profile,n_genes=len(x),covariate_rank_including_intercept=rank,
            unadjusted_spearman=float(stats.spearmanr(means[:,j],p0[:,j])[0]),
            partial_rank_correlation=float(np.corrcoef(residual.T)[0,1]),
            residual_rums_rank_sd=float(residual[:,0].std()), residual_probability_rank_sd=float(residual[:,1].std())))
    write_csv(args.out_dir/'multivariable_model_summary.csv', summary)
    write_csv(args.out_dir/'partial_rank_associations.csv',partial)
    write_csv(args.out_dir/'held_out_predictions.csv',prediction_rows)
    fig, axes = prior.paper.plt.subplots(1,2,figsize=(11,4))
    for name,_,_ in blocks:
        values=[next(r['gene_held_out_r_squared'] for r in summary if r['transform']=='log1p' and r['model']==name and r['Signature']==p) for p in profiles]
        axes[0].plot(profiles,values,'o-',label=name.replace('_',' '),markersize=4)
    axes[0].set(ylabel='Gene-held-out R-squared (log1p RUMS)')
    axes[0].legend(fontsize=7,frameon=False)
    axes[1].scatter(np.log1p(means[:,0]),log_predictions['plus_trinucleotide_opportunity'][:,0],s=3,alpha=.15)
    axes[1].set(xlabel='SBS2 observed log1p RUMS',ylabel='Held-out prediction from sequence covariates')
    prior.paper.save_plot(fig,args.out_dir,'composition_adjustment')
    report=dict(analysis_passed=True,stage='joint_composition_adjustment',fresh_simulation=False,
        representative_genes=len(x),profiles=profiles,primary_response='log1p RUMS',sensitivity_response='raw RUMS',
        covariates=names,models=[b[0] for b in blocks],folds=5,
        fit='Standardized ordinary least squares with SVD rank reduction; relative singular-value cutoff 1e-10',
        collinearity='GC/GC3, codon and context frequencies have structural dependencies; individual coefficients are not interpreted',
        validation='Fixed SHA256 gene folds shared across profiles; scaling and regression fitted only on training genes',
        partial_correlation='Pearson correlation between residual ranks of RUMS and weighted probability after projection on ranked covariates',
        causal_or_adaptive_inference=False,orthologue_or_species_adjustment=False,
        interpretation='Associations and predictive sequence-composition adjustment within the human representative-gene cohort, not a causal partition or a synonymous-recoding null.',
        input_manifest_sha256=digest(args.inputs/'input_manifest.json'),script_sha256=digest(__file__),
        focal_results=[r for r in summary if r['transform']=='log1p' and r['Signature']=='SBS2'],
        partial_rank_associations=partial,elapsed_seconds=time.perf_counter()-started,
        outputs=[dict(filename=p.name,sha256=digest(p)) for p in sorted(args.out_dir.iterdir())])
    (args.out_dir/'analysis_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('covariates','outputs')},indent=2),flush=True)


if __name__=='__main__':main()
