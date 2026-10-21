"""Explicit reconstructions of the original S11 sensitivity and S12 simulation graph."""
import argparse
import hashlib
import json
import multiprocessing as mp
import platform
from pathlib import Path
import re
import sys
import tempfile
import time

import numpy as np
import Bio
from Bio import SeqIO

from validate_manuscript_results import digest,require,write_csv,save_plot
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'02_run_sinabro'))
import run_rums as package_helpers


def initialize(source):
    package_helpers.load_package(Path(source))


def simulate(task):
    stage,identifier,sequence,first,count,seed,cap=task
    from sinabro import Trajectory
    rows=[];nodes={};edges=[]
    for index in range(first,first+count):
        key=hashlib.sha256(f'{seed}\0supplementary\0{identifier}\0{index}'.encode()).digest()
        path_seed=int.from_bytes(key[:8],'big')
        relaxed=Trajectory(identifier,sequence,seed=path_seed,cds_start=1,cds_end=len(sequence)-1)
        relaxed.autofill(condition='blosum',method='signature',mutational_signature='SBS2',
            cosmic_version='3.3.1',genome_ref='GRCh37',
            threshold=20,matrix_name='BLOSUM80',mode='global',gap_open=-10,gap_extend=-.5,
            stop_policy='score',max_steps=cap)
        require(relaxed.accepted_steps is not None,'Incomplete BLOSUM path; no mean will be published')
        row=dict(Transcript=identifier,path=index,seed=path_seed,blosum_neutral_steps=relaxed.accepted_steps,
                 terminal_reason=relaxed.terminal_reason,terminal_distance=relaxed.score)
        if stage=='s11':
            strict=Trajectory(identifier,sequence,seed=path_seed,cds_start=1,cds_end=len(sequence)-1)
            strict.autofill(condition='nonsynonymous',method='signature',mutational_signature='SBS2',
                            cosmic_version='3.3.1',genome_ref='GRCh37',max_steps=cap)
            require(strict.accepted_steps is not None,'Incomplete protein-identity path')
            require(relaxed.accepted_steps>=strict.accepted_steps,'Relaxed criterion stopped before strict identity')
            strict_records=strict.get_records();relaxed_records=relaxed.get_records()
            require(all(str(a.sequence)==str(b.sequence) for a,b in zip(strict_records,relaxed_records)),
                    'Matched criterion paths have different mutation prefixes')
            row['identity_neutral_steps']=strict.accepted_steps
        else:
            previous=None
            for step,record in enumerate(relaxed.get_records()):
                text=str(record.sequence)
                require(len(text)==len(sequence) and text[0]==sequence[0] and text[-1]==sequence[-1],
                        'Graph state changed sequence length or reference flanks')
                node=hashlib.sha256(text.encode()).hexdigest()
                differences=[f'{i}{old}>{new}' for i,(old,new) in enumerate(zip(sequence[1:-1],text[1:-1]),1) if old!=new]
                distance=float(record.note.split('score:',1)[1]) if record.note.startswith('score:') else 0.
                nodes[node]=dict(node=node,transcript_id=identifier,hamming_distance=len(differences),blosum_distance=distance,
                    cds_substitutions=';'.join(differences),observed_cosmic_status='not_available')
                if previous is not None:edges.append(dict(source=previous,target=node,path=index,step=step))
                previous=node
        rows.append(row)
    return rows,list(nodes.values()),edges


def audit_graph(rows,nodes,edges,reference):
    """Independently reconstruct graph DNA states and check each chronological path."""
    require(len(rows)==len({r['path'] for r in rows}),'Duplicate graph trajectory')
    root=hashlib.sha256(reference.encode()).hexdigest()
    require(root in nodes,'Graph lacks its reference state')
    substitutions={}
    for key,node in nodes.items():
        changes={};sequence=list(reference)
        for text in filter(None,node['cds_substitutions'].split(';')):
            match=re.fullmatch(r'([1-9][0-9]*)([ACGT])>([ACGT])',text)
            require(match is not None,'Invalid CDS substitution notation')
            position=int(match[1]);old,new=match[2],match[3]
            require(position<len(reference)-1 and position not in changes
                    and reference[position]==old and old!=new,'Invalid CDS substitution coordinate/reference')
            changes[position]=new;sequence[position]=new
        require(hashlib.sha256(''.join(sequence).encode()).hexdigest()==key==node['node'],
                'Graph state hash does not match its reference-relative substitutions')
        require(len(changes)==node['hamming_distance'],'Graph Hamming distance mismatch')
        require(node['observed_cosmic_status']=='not_available','Unverified observed-variant annotation')
        substitutions[key]=changes
    by_path={r['path']:[] for r in rows}
    for edge in edges:
        require(edge['path'] in by_path,'Unknown graph path')
        require(edge['source'] in nodes and edge['target'] in nodes,'Missing graph edge endpoint')
        left,right=substitutions[edge['source']],substitutions[edge['target']]
        difference=sum(left.get(p,reference[p])!=right.get(p,reference[p]) for p in set(left)|set(right))
        require(difference==1,'Graph edge is not a single-base substitution')
        by_path[edge['path']].append(edge)
    visited={root}
    for row in rows:
        path=sorted(by_path[row['path']],key=lambda e:e['step'])
        require(len(path)==row['blosum_neutral_steps']+1,'Graph path/event count mismatch')
        require([e['step'] for e in path]==list(range(1,len(path)+1)),
                'Missing or duplicate graph step')
        previous=root
        for edge in path:
            require(edge['source']==previous,'Disconnected graph trajectory')
            require(nodes[previous]['blosum_distance']<20,'Path continued beyond the BLOSUM threshold')
            previous=edge['target'];visited.add(previous)
        require(nodes[previous]['blosum_distance']==row['terminal_distance']>=20
                and row['terminal_reason']=='evaluation','Graph terminal state/criterion mismatch')
    require(visited==set(nodes),'Graph contains unused states')
    return dict(graph_audit_passed=True,complete_paths=len(rows),incomplete_paths=0,
                single_substitution_edges_checked=len(edges),reconstructed_dna_states_checked=len(nodes))


def run(args):
    require(args.workers>0 and args.n_traj>0 and args.max_steps>0 and args.seed>=0,'Invalid simulation settings')
    manifest=json.loads(args.cases.read_text())
    require(digest(args.input)==manifest['human_input_sha256'],'Supplementary FASTA differs from frozen input')
    chosen=manifest['s11'] if args.stage=='s11' else [manifest['s12']]
    if args.case_limit is not None:
        require(0<args.case_limit<=len(chosen),'Invalid bounded case limit')
        chosen=chosen[:args.case_limit]
    selected={r['transcript_id']:r for r in chosen}
    records={r.id.split('|')[0]:r for r in SeqIO.parse(args.input,'fasta') if r.id.split('|')[0] in selected}
    require(set(records)==set(selected),'Missing supplementary transcript')
    for key,record in records.items():
        require(record.id.split('|')[4]==selected[key]['transcript_name']
                and hashlib.sha256(str(record.seq).encode()).hexdigest()==selected[key]['sequence_sha256'],
                'Supplementary record changed')
    output=args.out_dir.resolve()
    require(not output.exists(),'Choose a new supplementary output directory')
    require(args.input.resolve().parent not in output.parents,'Output nested inside input directory')
    package_helpers.load_package(args.sinabro_source)
    source=package_helpers.source
    package_sha=package_helpers.package_digest()
    signature_file=package_helpers.signature_path('3.3.1','GRCh37')
    started=time.perf_counter()
    report=dict(stage=args.stage,simulation_passed=False,
        historical_subset_recovered=(args.stage=='s12' and chosen[0].get('historical_transcript_recovered',False)),
        original_figure_reproduced=False,observed_cosmic_variant_overlay_reproduced=False,
        fresh_simulation=True,corrected_package_blosum_recalculation=True,
        input_sha256=digest(args.input),cases_sha256=digest(args.cases),package_source_sha256=package_sha,
        signature_version='3.3.1',signature_reference='GRCh37',signature_sha256=digest(signature_file),
        script_sha256=digest(__file__),signatures=['SBS2'],n_trajectories=args.n_traj,seed=args.seed,
        selected_transcripts=list(selected),case_limit=args.case_limit,max_steps=args.max_steps,
        bounded_software_check=args.case_limit is not None,
        runtime=dict(python=platform.python_version(),numpy=np.__version__,biopython=Bio.__version__),
        criterion=dict(matrix='BLOSUM80',mode='global',gap_open=-10,gap_extend=-.5,threshold=20,
            threshold_rule='distance >= 20',stop_policy='score'),
        interpretation='Fresh corrected-package simulation of the documented cases. Historical random streams and observed COSMIC annotations were not recovered. BLOSUM distance is a proxy, not experimentally measured functional loss.')
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.paper-supplement-',dir=output.parent) as temporary:
        stage=Path(temporary)/'results';stage.mkdir()
        tasks=[(args.stage,key,str(records[key].seq),first,min(25,args.n_traj-first),args.seed,args.max_steps)
               for key in selected for first in range(0,args.n_traj,25)]
        values=[];nodes={};edges=[]
        with mp.get_context('spawn').Pool(args.workers,initializer=initialize,initargs=(str(source),)) as pool:
            for i,(rows,vertices,connections) in enumerate(pool.imap_unordered(simulate,tasks),1):
                values.extend(rows);edges.extend(connections)
                for vertex in vertices:
                    if vertex['node'] in nodes:
                        require(vertex==nodes[vertex['node']],'Inconsistent graph state score')
                    nodes[vertex['node']]=vertex
                if i%10==0 or i==len(tasks):print(f'{args.stage}: finished {i}/{len(tasks)} path blocks',flush=True)
        values.sort(key=lambda row:(row['Transcript'],row['path']))
        require(len(values)==len(selected)*args.n_traj,'Wrong supplementary path count')
        write_csv(stage/'trajectory_counts.csv',values)
        if args.stage=='s11':
            summary=[]
            for key,annotation in selected.items():
                rows=[r for r in values if r['Transcript']==key]
                strict=np.array([r['identity_neutral_steps'] for r in rows]);relaxed=np.array([r['blosum_neutral_steps'] for r in rows])
                summary.append(dict(**annotation,n_trajectories=len(rows),identity_mean=float(strict.mean()),
                    blosum_mean=float(relaxed.mean()),paired_difference_mean=float((relaxed-strict).mean()),
                    paired_difference_mc_se=float((relaxed-strict).std(ddof=1)/np.sqrt(len(rows))) if len(rows)>1 else None))
            write_csv(stage/'s11_criterion_comparison.csv',summary)
            fig,ax=plt.subplots(figsize=(6,4))
            ax.scatter([r['identity_mean'] for r in summary],[r['blosum_mean'] for r in summary])
            limit=max(max(r['identity_mean'],r['blosum_mean']) for r in summary)
            ax.plot([0,limit],[0,limit],'--',color='grey')
            ax.set(xlabel='Protein-identity neutral steps',ylabel='BLOSUM80 criterion-neutral steps',title='S11 reconstruction: frozen gene subset')
            save_plot(fig,stage,'figureS11_reconstructed')
            report['summary']=summary
            report['matched_path_prefixes_checked']=True
        else:
            reference=str(records[next(iter(selected))].seq)
            report.update(audit_graph(values,nodes,edges,reference))
            record=chosen[0]
            (stage/'s12_reference.cdsplus.fa').write_text(
                f'>{record["transcript_id"]}|{record["transcript_name"]}|cds_start_0based=1|cds_end_exclusive={len(reference)-1}\n{reference}\n')
            ordered=sorted(nodes)
            coordinates={key:(0.,0.) if nodes[key]['hamming_distance']==0 else (
                nodes[key]['hamming_distance']*np.cos(2*np.pi*int(key[:8],16)/2**32),
                nodes[key]['hamming_distance']*np.sin(2*np.pi*int(key[:8],16)/2**32)) for key in ordered}
            write_csv(stage/'s12_graph_nodes.csv',[dict(nodes[key],x=coordinates[key][0],y=coordinates[key][1]) for key in ordered])
            write_csv(stage/'s12_graph_edges.csv',sorted(edges,key=lambda row:(row['path'],row['step'])))
            fig,ax=plt.subplots(figsize=(6,5))
            ax.add_collection(LineCollection([[coordinates[e['source']],coordinates[e['target']]] for e in edges],alpha=.2,linewidth=.4))
            scatter=ax.scatter([coordinates[k][0] for k in ordered],[coordinates[k][1] for k in ordered],
                c=[nodes[k]['blosum_distance'] for k in ordered],s=8,cmap='viridis')
            fig.colorbar(scatter,ax=ax,label='BLOSUM80 alignment-score distance')
            ax.set_aspect('equal');ax.autoscale_view()
            ax.set_title('BRCA1-203 simulation graph\nObserved COSMIC overlay unavailable',fontsize=10)
            ax.set(xlabel='Radial distance = DNA Hamming distance',ylabel='Deterministic angular layout')
            save_plot(fig,stage,'figureS12_simulation_only')
            report.update(nodes=len(nodes),edges=len(edges),
                reference_sequence_sha256=hashlib.sha256(reference.encode()).hexdigest(),
                coordinate_system='CDS-relative 1-based substitutions against the initial reference; flanking bases excluded',
                graph_layout='Deterministic angle, radius equals reference DNA Hamming distance',
                historical_script_n_sim=record.get('historical_script_n_sim'),
                blosum_neutral_steps_mean=float(np.mean([r['blosum_neutral_steps'] for r in values])),
                blosum_neutral_steps_median=float(np.median([r['blosum_neutral_steps'] for r in values])),
                missing_input='Observed COSMIC breast-cancer variants with matched BRCA1 transcript coordinates')
        require(digest(args.input)==report['input_sha256'] and package_helpers.package_digest()==package_sha
                and digest(args.cases)==report['cases_sha256']
                and digest(signature_file)==report['signature_sha256'],
                'Scientific inputs changed during supplementary run')
        report.update(simulation_passed=True,elapsed_seconds=time.perf_counter()-started,
            outputs=[dict(filename=p.name,sha256=digest(p)) for p in sorted(stage.iterdir())])
        (stage/'supplementary_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        require(not output.exists(),'Supplementary output appeared during computation')
        stage.rename(output)
    print(json.dumps({k:v for k,v in report.items() if k!='outputs'},indent=2),flush=True)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('s11','s12'),required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--sinabro-source',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--cases',type=Path,default=HERE/'supplementary_cases.json')
    parser.add_argument('--n-traj',type=int,default=1000)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--max-steps',type=int,default=10000)
    parser.add_argument('--case-limit',type=int,help='Explicit bounded software smoke test')
    run(parser.parse_args(argv))


if __name__=='__main__':
    main()
