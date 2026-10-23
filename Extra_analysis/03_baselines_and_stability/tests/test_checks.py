"""Brute-force sequence oracles and independent saved-path numerical checks."""
from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from Bio.Seq import Seq

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import analyze_baselines as baseline
import analyze_stability as stability
from common import prior, digest, write_csv


def channels():
    return [f'{a}[{ref}>{alt}]{b}' for ref in 'CT' for alt in 'ACGT' if alt!=ref
            for a in 'ACGT' for b in 'ACGT']


class NumericalTests(unittest.TestCase):
    def test_all_codon_uniform_counts_match_full_translation_enumeration(self):
        lookup=baseline.channel_lookup(channels())
        for codon in baseline.CODONS:
            sequence='C'+codon+'G'
            total,neutral=baseline.opportunities(sequence,lookup)
            original=str(Seq(codon).translate())
            count=0
            for i in range(3):
                for alt in 'ACGT':
                    if alt!=codon[i]:
                        mutant=codon[:i]+alt+codon[i+1:]
                        count+=str(Seq(mutant).translate())==original
            self.assertEqual(total.sum(),9)
            self.assertEqual(neutral.sum(),count)

    def test_weighted_channel_counts_match_package_event_probabilities(self):
        source=HERE.parents[2]/'sinabro'
        if not (source/'sinabro/core.py').exists():
            self.skipTest('Independent Sinabro package oracle requires the sibling checkout')
        sys.path.insert(0,str(source))
        from sinabro.core import probability_matrix,full_signature
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'signature.tsv'
            labels=channels();weights=np.random.default_rng(7).random(96)
            path.write_text('Type\ttest\n'+''.join(f'{s}\t{v}\n' for s,v in zip(labels,weights)))
            lookup=baseline.channel_lookup(labels)
            for seq in ['CATGAAATAAG','GATGGGCTAAC','CTGATAGTAAG','AAAAAAAT']:
                counts,neutral=baseline.opportunities(seq,lookup)
                actual=neutral@weights/(counts@weights)
                p=probability_matrix(seq,full_signature(path,'test'))
                protein=str(Seq(seq[1:-1]).translate());expected=0.
                for i in range(1,len(seq)-1):
                    for j,alt in enumerate('ACGT'):
                        if alt!=seq[i]:
                            mutant=seq[:i]+alt+seq[i+1:]
                            if str(Seq(mutant[1:-1]).translate())==protein:expected+=p[j,i]
                self.assertAlmostEqual(actual,expected,places=13)

    def test_bad_frame_and_reverse_complement_pair_counts(self):
        lookup=baseline.channel_lookup(channels())
        np.testing.assert_array_equal(np.bincount(lookup[lookup>=0]),np.full(96,2))
        for text in ['CATGAG','CATNNNTAGG']:
            with self.assertRaises(ValueError):baseline.opportunities(text,lookup)

    def test_prefix_metrics_match_direct_formulas_including_zero_paths(self):
        counts=np.array([[0,2],[1,0],[2,3],[1,1]],dtype=np.int32)
        result=stability.prefix_statistics(counts,[2,4])
        np.testing.assert_allclose(result[0]['mean'],[.5,1])
        np.testing.assert_allclose(result[0]['zero_fraction'],[.5,.5])
        np.testing.assert_allclose(result[0]['mc_se'],counts[:2].std(axis=0,ddof=1)/np.sqrt(2))
        np.testing.assert_allclose(result[0]['delta'],counts[:2].mean(axis=0)-counts.mean(axis=0))
        np.testing.assert_array_equal(result[-1]['delta'],[0,0])
        with self.assertRaises(ValueError):stability.prefix_statistics(counts,[2,3])


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.export=self.root/'export';self.export.mkdir()
        self.runs=self.root/'run';(self.runs/'transcripts').mkdir(parents=True)
        self.audit=self.root/'audit';self.audit.mkdir()
        (self.root/'sequences').mkdir();self.fasta=self.root/'sequences/input.fa'
        (self.root/'sources').mkdir();self.signatures=self.root/'sources/signatures.tsv'
        names=['SBS1','SBS2','SBS13']
        labels=channels();weights=np.random.default_rng(1).random((96,3))
        self.signatures.write_text('Type\t'+'\t'.join(names)+'\n'+''.join(s+'\t'+'\t'.join(map(str,v))+'\n' for s,v in zip(labels,weights)))
        annotation=[];seqs=[]
        for i in range(24):
            identifier=f'ENST{i:05d}.1';gene=f'GENE{i//2}.1'
            cds='ATG'+['AAA','GGC','TTC'][i%3]*(i+1)+'TAA'
            fid=f'{identifier}|{gene}|x|y|name{i}|g{i//2}|'
            seqs.append(f'>{fid}\nC{cds}G\n')
            annotation.append(dict(transcript_id=identifier,gene_id=gene,gene_name=f'g{i//2}',
                cds_nt=len(cds),fasta_id=fid,internal_stop='False'))
        self.fasta.write_text(''.join(seqs))
        config=dict(runner='rums-v3',translation_table=1,cds_start=1,cds_end='len(sequence)-1',sinabro='3.0.0',
            criterion='protein_identity',n_trajectories=1000,signatures=names,max_steps=10000,master_seed=42,
            input_sha256=digest(self.fasta),signature_sha256=digest(self.signatures))
        fingerprint=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
        (self.runs/'run_manifest.json').write_text(json.dumps(dict(configuration=config,config_sha256=fingerprint)))
        rows=[];means=[]
        rng=np.random.default_rng(42)
        for i,a in enumerate(annotation):
            counts=rng.poisson(np.array([.4,.8,.2])*(1+i/10),size=(1000,3)).astype(np.int32)
            metadata=dict(config_sha256=fingerprint,transcript_id=a['transcript_id'],fasta_id=a['fasta_id'],
                          signatures=names,n_trajectories=1000)
            path=self.runs/'transcripts'/(a['transcript_id']+'.npz')
            np.savez_compressed(path,accepted_steps=counts,terminal_codes=np.zeros_like(counts,dtype=np.uint8),
                n_complete=np.full(3,1000),n_incomplete=np.zeros(3,dtype=int),metadata=json.dumps(metadata))
            rows.append(dict(transcript_id=a['transcript_id'],sha256=digest(path)))
            means.append(dict(Transcript=a['transcript_id'],**dict(zip(names,counts.mean(axis=0)))))
        write_csv(self.audit/'transcript_audit.csv',rows)
        audit=dict(validation_passed=True,all_paths_complete=True,issues=[],configuration=config,config_sha256=fingerprint,
                   manifest_sha256=digest(self.runs/'run_manifest.json'))
        (self.audit/'validation_summary.json').write_text(json.dumps(audit))
        write_csv(self.export/'rums_mean_by_transcript_id.csv',means)
        write_csv(self.export/'transcript_annotations.csv',annotation)
        manifest=dict(aggregation_passed=True,configuration=config,config_sha256=fingerprint,transcripts=24,
            checkpoint_hashes_matched=24,input_path=str(self.fasta),input_sha256=digest(self.fasta),
            run_directory=str(self.runs),validation_directory=str(self.audit),
            validation_report_sha256=digest(self.audit/'validation_summary.json'),
            audit_inventory_sha256=digest(self.audit/'transcript_audit.csv'),
            outputs=[dict(filename=p.name,sha256=digest(p)) for p in sorted(self.export.iterdir())])
        (self.export/'aggregation_manifest.json').write_text(json.dumps(manifest))

    def test_both_full_prefix_workflows_leave_inputs_unchanged(self):
        before={p:digest(p) for p in self.root.rglob('*') if p.is_file()}
        args=['--human-export',str(self.export),'--allow-subset']
        with redirect_stdout(io.StringIO()):
            a=baseline.main([*args,'--signatures',str(self.signatures),'--out-dir',str(self.root/'baseline')])
            b=stability.main([*args,'--per-bin','2','--out-dir',str(self.root/'stability')])
        self.assertTrue(a['analysis_passed']);self.assertTrue(b['analysis_passed'])
        self.assertFalse(a['fresh_simulation']);self.assertFalse(b['fresh_simulation'])
        self.assertEqual(a['representative_genes'],12)
        self.assertEqual(b['prefixes'],[100,200,500,1000])
        self.assertEqual(b['prefix_summary'][-1]['maximum_absolute_delta_from_full'],0)
        self.assertFalse(b['universal_1000_path_sufficiency_claimed'])
        self.assertTrue(all(digest(p)==sha for p,sha in before.items()))
        self.assertTrue((self.root/'baseline/one_step_baseline_comparison.pdf').read_bytes().startswith(b'%PDF'))

    def test_sampling_is_deterministic_and_independent_of_rums(self):
        data=prior.load_dataset(self.export);reps=prior.representatives(data)
        selected,_=stability.select_sample(data,reps,1,42)
        data['matrix'][:]=np.random.default_rng(99).normal(size=data['matrix'].shape)
        changed,_=stability.select_sample(data,reps,1,42)
        self.assertEqual(selected,changed)
        self.assertEqual(len({r['gene_id'] for r in selected}),len(selected))

    def test_changed_checkpoint_or_audit_is_rejected(self):
        path=next((self.runs/'transcripts').glob('*.npz'))
        data=prior.load_dataset(self.export)
        i=data['ids'].index(path.stem);expected=digest(path)
        path.write_bytes(path.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'changed since full audit'):
            stability.load_paths(path,expected,data,i)
        p=self.audit/'transcript_audit.csv';p.write_text(p.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'provenance mismatch'):
            stability.main(['--human-export',str(self.export),'--allow-subset','--out-dir',str(self.root/'bad')])


if __name__=='__main__':unittest.main()
