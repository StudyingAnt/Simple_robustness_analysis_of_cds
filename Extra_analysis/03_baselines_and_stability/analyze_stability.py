"""Reuse saved trajectories to describe Monte Carlo precision and prefix stability."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import time

from common import (np, stats, prior, require, digest, write_csv, correlation,
                    load_human, check_output, finish, add_arguments)


def select_sample(data, representatives, per_bin, seed):
    require(per_bin > 0 and seed >= 0, 'Invalid sampling settings')
    length = np.array([data['features'][i]['cds_nt'] for i in representatives])
    gc3 = np.array([data['features'][i]['gc3'] for i in representatives])
    cut_length, cut_gc3 = np.quantile(length, [1/3,2/3]), np.quantile(gc3, [1/3,2/3])
    bins = {}
    for n,i in enumerate(representatives):
        key = (int(np.searchsorted(cut_length, length[n], side='right')),
               int(np.searchsorted(cut_gc3, gc3[n], side='right')))
        bins.setdefault(key, []).append(i)
    selected, strata = set(), {}
    for key, indices in sorted(bins.items()):
        ordered = sorted(indices, key=lambda i: hashlib.sha256(f'{seed}:{data["ids"][i]}'.encode()).hexdigest())
        selected.update(ordered[:per_bin])
        strata[key] = dict(available=len(indices), sampled=min(per_bin,len(indices)))
    # Include extremes by metadata alone, without inspecting RUMS or path outcomes.
    extremes = set()
    for array in (length,gc3):
        extremes.update([representatives[int(np.argmin(array))],representatives[int(np.argmax(array))]])
    selected.update(extremes)
    rows = []
    for i in sorted(selected, key=lambda i:data['ids'][i]):
        f = data['features'][i]
        key = (int(np.searchsorted(cut_length,f['cds_nt'],side='right')),
               int(np.searchsorted(cut_gc3,f['gc3'],side='right')))
        rows.append(dict(index=i,Transcript=data['ids'][i],gene_id=data['annotations'][i]['gene_id'],
                         cds_nt=f['cds_nt'],gc3=f['gc3'],length_tertile=key[0]+1,gc3_tertile=key[1]+1,
                         metadata_extreme=i in extremes))
    return rows, dict(length_cutpoints=list(map(float,cut_length)),gc3_cutpoints=list(map(float,cut_gc3)),
                      strata={str(k):v for k,v in strata.items()})


def load_paths(path, expected_sha, data, index):
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected_sha, 'Checkpoint changed since full audit')
    with np.load(io.BytesIO(raw), allow_pickle=False) as z:
        meta = json.loads(str(z['metadata'].item()))
        require(meta['config_sha256'] == data['manifest']['config_sha256']
                and meta['transcript_id'] == data['ids'][index]
                and meta['fasta_id'] == data['annotations'][index]['fasta_id']
                and meta['signatures'] == data['signatures']
                and meta['n_trajectories'] == data['manifest']['configuration']['n_trajectories'],
                'Checkpoint metadata differs from export')
        counts, codes = z['accepted_steps'], z['terminal_codes']
        require(counts.shape == codes.shape == (meta['n_trajectories'],len(data['signatures']))
                and counts.dtype == np.int32 and codes.dtype == np.uint8
                and (codes == 0).all() and (counts >= 0).all()
                and (counts < data['manifest']['configuration']['max_steps']).all(),
                'Invalid or incomplete saved trajectories')
        require(np.array_equal(z['n_complete'],np.full(len(data['signatures']),meta['n_trajectories']))
                and (z['n_incomplete'] == 0).all()
                and np.allclose(counts.mean(axis=0),data['matrix'][index],rtol=0,atol=1e-12),
                'Checkpoint counts do not reproduce exported means')
    return counts


def prefix_statistics(counts, prefixes):
    require(prefixes == sorted(set(prefixes)) and prefixes[-1] == len(counts)
            and min(prefixes) >= 2 and max(prefixes) <= len(counts), 'Invalid prefix sizes')
    result = []
    reference = counts.mean(axis=0)
    for n in prefixes:
        values = counts[:n]
        means = values.mean(axis=0)
        sd = values.std(axis=0,ddof=1)
        result.append(dict(n=n,mean=means,sd=sd,mc_se=sd/np.sqrt(n),
                           zero_fraction=(values==0).mean(axis=0),delta=means-reference))
    return result


def analyze(args):
    started = time.perf_counter()
    data, representatives = load_human(args)
    prefixes = sorted(set(args.prefixes))
    require(len(prefixes) == len(args.prefixes) and min(prefixes)>=2 and prefixes[-1] == args.n_traj,
            'Prefixes must be distinct, >=2, and end at the saved trajectory count')
    manifest = data['manifest']
    runs = Path(args.run_dir or manifest['run_directory']).resolve()
    audit = Path(args.audit_dir or manifest['validation_directory']).resolve()
    report_file, inventory_file = audit/'validation_summary.json', audit/'transcript_audit.csv'
    report_hash, inventory_hash = digest(report_file), digest(inventory_file)
    require(report_hash == manifest['validation_report_sha256']
            and inventory_hash == manifest['audit_inventory_sha256'], 'Audit/export provenance mismatch')
    audit_report = json.loads(report_file.read_text())
    run_manifest = json.loads((runs/'run_manifest.json').read_text())
    require(audit_report['validation_passed'] is True and audit_report['all_paths_complete'] is True
            and audit_report['issues'] == [] and audit_report['configuration'] == manifest['configuration']
            and audit_report['config_sha256'] == manifest['config_sha256']
            and run_manifest['configuration'] == manifest['configuration']
            and run_manifest['config_sha256'] == manifest['config_sha256']
            and digest(runs/'run_manifest.json') == audit_report['manifest_sha256'], 'Unverified run configuration')
    with inventory_file.open(newline='') as stream:
        inventory = list(csv.DictReader(stream))
    hashes = {r['transcript_id']:r['sha256'] for r in inventory}
    require(len(inventory) == len(hashes) == len(data['ids']) and set(hashes) == set(data['ids']),
            'Audit inventory does not match export cohort')
    selected, strata = select_sample(data,representatives,args.per_bin,args.sample_seed)
    output = check_output(args.out_dir,[data['folder'],data['fasta'].parent,runs,audit])
    paths_by_prefix = {n:[] for n in prefixes}
    with tempfile.TemporaryDirectory(prefix='.path-stability-',dir=output.parent) as temp:
        stage = Path(temp)/'results';stage.mkdir()
        write_csv(stage/'sampled_representatives.csv',[{k:v for k,v in r.items() if k!='index'} for r in selected])
        detailed=[]
        for completed,row in enumerate(selected,1):
            i=row['index'];identifier=row['Transcript']
            counts=load_paths(runs/'transcripts'/(identifier+'.npz'),hashes[identifier],data,i)
            for values in prefix_statistics(counts,prefixes):
                paths_by_prefix[values['n']].append(values)
                for j,name in enumerate(data['signatures']):
                    detailed.append(dict(Transcript=identifier,Signature=name,n_paths=values['n'],
                        mean_rums=float(values['mean'][j]),path_sd=float(values['sd'][j]),mc_se=float(values['mc_se'][j]),
                        zero_path_fraction=float(values['zero_fraction'][j]),
                        mean_minus_full_saved_mean=float(values['delta'][j])))
            if completed%25==0 or completed==len(selected):
                print(f'Checked saved paths {completed}/{len(selected)} transcripts',flush=True)
        write_csv(stage/'path_precision_by_case_and_prefix.csv',detailed)
        summary=[];profile_rows=[]
        full_means=np.array([v['mean'] for v in paths_by_prefix[args.n_traj]])
        reference_profile_mean=full_means.mean(axis=0)
        for n in prefixes:
            values=paths_by_prefix[n]
            sem=np.array([v['mc_se'] for v in values]);delta=np.array([v['delta'] for v in values])
            zero=np.array([v['zero_fraction'] for v in values])
            means=np.array([v['mean'] for v in values]);profile_mean=means.mean(axis=0)
            summary.append(dict(n_paths=n,cases=sem.size,median_mc_se=float(np.median(sem)),
                p90_mc_se=float(np.quantile(sem,.9)),maximum_mc_se=float(sem.max()),
                median_absolute_delta_from_full=float(np.median(np.abs(delta))),
                p90_absolute_delta_from_full=float(np.quantile(np.abs(delta),.9)),
                maximum_absolute_delta_from_full=float(np.abs(delta).max()),
                mean_zero_path_fraction=float(zero.mean()),
                profile_mean_rank_spearman_vs_full=correlation(profile_mean,reference_profile_mean),
                sample_highest_mean_profile=data['signatures'][int(np.argmax(profile_mean))]))
            ranks=stats.rankdata(-profile_mean,method='min')
            for j,name in enumerate(data['signatures']):
                profile_rows.append(dict(n_paths=n,Signature=name,sampled_transcripts=len(selected),
                    mean_across_sample=float(profile_mean[j]),rank_mean_desc=float(ranks[j]),
                    median_mc_se=float(np.median(sem[:,j])),p90_mc_se=float(np.quantile(sem[:,j],.9)),
                    mean_zero_path_fraction=float(zero[:,j].mean()),
                    p90_absolute_delta_from_full=float(np.quantile(np.abs(delta[:,j]),.9))))
        write_csv(stage/'precision_summary_by_prefix.csv',summary)
        write_csv(stage/'precision_by_signature_and_prefix.csv',profile_rows)
        fig,axes=prior.paper.plt.subplots(1,2,figsize=(9,4))
        axes[0].plot(prefixes,[r['median_mc_se'] for r in summary],'o-',label='Median')
        axes[0].plot(prefixes,[r['p90_mc_se'] for r in summary],'o-',label='90th percentile')
        axes[0].set(xlabel='Number of saved paths used',ylabel='Empirical Monte Carlo SE');axes[0].legend()
        axes[1].plot(prefixes,[r['p90_absolute_delta_from_full'] for r in summary],'o-')
        axes[1].set(xlabel='Number of saved paths used',ylabel='90th percentile absolute difference from full saved mean')
        prior.paper.save_plot(fig,stage,'saved_path_precision')
        require(digest(report_file)==report_hash and digest(inventory_file)==inventory_hash
                and digest(data['fasta'])==manifest['input_sha256'], 'Input changed during analysis')
        report=dict(analysis_passed=True,stage='saved_path_stability',fresh_simulation=False,
            bounded_export=args.allow_subset,stability_sample_is_full_genome=False,
            saved_paths_per_case=args.n_traj,representative_pool=len(representatives),
            sampled_transcripts=len(selected),signatures=len(data['signatures']),prefixes=prefixes,
            sample_seed=args.sample_seed,per_feature_bin=args.per_bin,stratification=strata,
            selection='Length/GC3 tertiles, SHA256 ordering within each bin; metadata extrema included; no RUMS-based selection',
            export_manifest_sha256=data['manifest_sha256'],audit_report_sha256=report_hash,
            audit_inventory_sha256=inventory_hash,selected_checkpoints_sha256={r['Transcript']:hashes[r['Transcript']] for r in selected},
            script_sha256=digest(__file__),common_sha256=digest(Path(__file__).parent/'common.py'),
            prefix_summary=summary,focal_profiles=[r for r in profile_rows if r['Signature'] in ('SBS2','SBS13','SBS19')],
            interpretation='Nested prefixes describe sensitivity within existing draws; the full saved mean is not independent ground truth. Empirical MC SE is distinct from gene/species variation. Zero empirical variance does not establish zero true variance.',
            universal_1000_path_sufficiency_claimed=False,independent_convergence_to_truth_test=False,
            elapsed_seconds=time.perf_counter()-started)
        finish(stage,output,report)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    parser.add_argument('--run-dir',type=Path,help='Optional relocated run with identical manifest and checkpoint hashes')
    parser.add_argument('--audit-dir',type=Path)
    parser.add_argument('--prefixes',type=int,nargs='+',default=[100,200,500,1000])
    parser.add_argument('--per-bin',type=int,default=20)
    parser.add_argument('--sample-seed',type=int,default=42)
    return analyze(parser.parse_args(argv))


if __name__=='__main__':
    main()
