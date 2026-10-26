"""Independent numerical checks and bounded original-paper pipeline exercise."""
import argparse
from contextlib import redirect_stdout
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy import stats

PAPER = Path(__file__).resolve().parents[1] / '07_robustness_analysis'
sys.path.insert(0, str(PAPER))
import analyze_robustness_results as paper
import validate_signature_sampler as sampler


def channels():
    return [a + '['+ref+'>'+alt+']'+b for ref in 'CT' for alt in 'ACGT' if alt != ref
            for a in 'ACGT' for b in 'ACGT']


class NumericalTests(unittest.TestCase):
    def test_pca_reconstructs_centered_matrix_and_covariance_eigenvalues(self):
        x = np.random.default_rng(3).normal(size=(30, 5)) * np.arange(1, 6)
        scores, components, ratios = paper.raw_pca(x)
        np.testing.assert_allclose(scores @ components, x-x.mean(axis=0), atol=1e-12)
        eigenvalues = np.linalg.eigvalsh(np.cov(x, rowvar=False))[::-1]
        np.testing.assert_allclose(ratios, eigenvalues / eigenvalues.sum(), atol=1e-12)
        self.assertTrue(np.all(np.diff(ratios) <= 0))

    def test_dunn_two_groups_matches_independent_mannwhitney_without_continuity(self):
        groups = [np.array([1,2,2,3]), np.array([2,4,5,6])]
        result = paper.dunn(groups, ['a','b'])[0]
        reference = stats.mannwhitneyu(*groups, method='asymptotic', use_continuity=False)
        self.assertAlmostEqual(result['p_two_sided'], reference.pvalue, places=13)
        result = paper.dunn(groups + [np.array([1,5,8])], ['a','b','c'])
        self.assertEqual(len(result), 3)
        self.assertTrue(all(r['p_bonferroni'] == min(1.,3*r['p_two_sided']) for r in result))

    def test_zero_probabilities_have_defined_entropy_and_feature_fractions(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'signatures.tsv'
            labels = channels()
            with path.open('w') as stream:
                stream.write('Type\tuniform\tpoint\n')
                for i,label in enumerate(labels):stream.write(f'{label}\t{1/96}\t{int(i==0)}\n')
            features = paper.signature_features(path,['uniform','point'])
            self.assertAlmostEqual(features['uniform']['entropy_bits'],np.log2(96))
            self.assertAlmostEqual(features['uniform']['transition_fraction'],1/3)
            self.assertAlmostEqual(features['uniform']['gc_target_preference'],.5)
            self.assertEqual(features['point']['entropy_bits'],0)

    def test_gc3_uses_coding_frame_and_excludes_flanks(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'input.fa'
            path.write_text('>id\nCATGGGCTAAG\n')
            row = paper.gc_features(path,['id'])[0]
            self.assertEqual(row['cds_nt'],9)
            self.assertAlmostEqual(row['gc3'],2/3)
            self.assertAlmostEqual(row['gc_total'],4/9)

    def test_conditional_sampler_oracle_respects_sequence_opportunities(self):
        labels = channels()
        probabilities = sampler.conditional_channel_probabilities('AAAAAA',labels,np.ones(96)/96)
        nonzero = [labels[i] for i in np.flatnonzero(probabilities)]
        self.assertEqual(set(nonzero),{'T[T>A]T','T[T>C]T','T[T>G]T'})
        np.testing.assert_allclose(probabilities[probabilities>0],np.full(3,1/3))

    def test_sampler_rejects_wrong_profile_even_in_rare_channels(self):
        checks,_ = sampler.channel_checks(np.array([100,0]),np.array([.001,99.999]),
                                          np.array([.001,.001]),100)
        self.assertFalse(checks['passed'])
        json.dumps(checks,allow_nan=False)

    @unittest.skipUnless((PAPER.parent.parent/'sinabro/sinabro/core.py').exists(), 'Requires sibling package source')
    def test_sampler_spawn_cli_uses_requested_profile_and_writes_serializable_report(self):
        with tempfile.TemporaryDirectory() as temp:
            output=Path(temp)/'output'
            result=subprocess.run([sys.executable,str(PAPER/'validate_signature_sampler.py'),
                '--sinabro-source',str(PAPER.parent.parent/'sinabro'),'--out-dir',str(output),
                '--n-sequences','40','--length','100','--workers','2','--signatures','SBS2','SBS13'],
                capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            report=json.loads((output/'sampler_validation_summary.json').read_text())
            self.assertTrue(report['validation_passed'])
            self.assertEqual(report['total_single_substitutions'],80)
            self.assertTrue(all(row['cosine_to_original_profile']>.7 for row in report['results']))


class PipelineTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.signature = self.root/'signatures.tsv'
        with self.signature.open('w') as stream:
            stream.write('Type\tSBS1\tSBS2\tSBS13\n')
            labels=channels()
            weights=np.array([[1, 1+10*(label[2:5] in ('C>T','T>C')),
                               1+20*(label[2:5] not in ('C>T','T>C'))] for label in labels],dtype=float)
            weights/=weights.sum(axis=0)
            for label,row in zip(labels,weights):stream.write(label+'\t'+'\t'.join(map(str,row))+'\n')
        self.source_hash = paper.digest(self.signature)
        self.cosmic = ['SBS1','SBS2','SBS13']
        self.enzyme = [f'HL{s}_{e}' for s in ('molMol2','myoMyo6') for e in
                      (('APOBEC3A','APOBEC3C','APOBEC1') if s=='molMol2' else ('APOBEC3A3','APOBEC3C','APOBEC1'))]
        self.fasta = {}
        self.matrices = {}
        for dataset in paper.SPECIES:
            folder = self.root/'sequences'/dataset
            folder.mkdir(parents=True)
            path = folder/'input.fa'
            path.write_text(''.join(f'>{dataset}{i}\nCATG'+('AAA','GGC','CTG','TTC','GAA')[i]*(i+1)+'TAAG\n' for i in range(5)))
            self.fasta[dataset] = path
            self.make_export(self.root/'apobec'/dataset,dataset,self.enzyme)
            self.make_export(self.root/'human' if dataset=='human' else self.root/'bats'/dataset,dataset,self.cosmic)
        self.args=argparse.Namespace(human_export=self.root/'human',bat_export=self.root/'bats',
            apobec_export=self.root/'apobec',signatures=self.signature,out_dir=self.root/'results',
            n_traj=10,allow_subset=True)

    def make_export(self,folder,dataset,names):
        folder.mkdir(parents=True)
        rng=np.random.default_rng(len(names)+paper.SPECIES.index(dataset))
        matrix=rng.integers(0,100,size=(5,len(names)))/10
        self.matrices[str(folder)]=matrix
        path=folder/'rums_mean_by_transcript_id.csv'
        with path.open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['Transcript',*names])
            writer.writerows([[dataset+str(i),*row] for i,row in enumerate(matrix)])
        config=dict(signatures=names,criterion='protein_identity',n_trajectories=10,
                    signature_sha256=self.source_hash,input_sha256=paper.digest(self.fasta[dataset]),
                    mutation_model='experimental_apobec_NNC_v1' if names==self.enzyme else 'COSMIC')
        manifest=dict(aggregation_passed=True,configuration=config,
            config_sha256=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest(),
            transcripts=5,checkpoint_hashes_matched=5,input_path=str(self.fasta[dataset]),
            input_sha256=config['input_sha256'],outputs=[dict(filename=path.name,sha256=paper.digest(path))])
        (folder/'aggregation_manifest.json').write_text(json.dumps(manifest))

    def test_bounded_end_to_end_preserves_direction_and_writes_all_main_panels(self):
        before={str(p):paper.digest(p) for p in self.root.rglob('*') if p.is_file()}
        with redirect_stdout(io.StringIO()):report=paper.analyze(self.args)
        self.assertTrue(report['analysis_passed'])
        self.assertTrue(report['bounded_subset'])
        self.assertEqual(set(report['pca']),{'all','noseqerr','aetiology'})
        self.assertEqual(len(report['native_panel_tests']),3)
        self.assertEqual(report['native_panel_tests'][1]['original_direction'],'less')
        for name in ('figure3_pca','figure4_pc_loadings','figure5_signature_features',
                     'figure7g_native_processes','figureS9_experimental_enzymes'):
            self.assertTrue((self.args.out_dir/(name+'.pdf')).read_bytes().startswith(b'%PDF'))
        expected = self.matrices[str(self.root/'human')][:,1].mean()
        self.assertEqual(report['native_panel']['human']['mean'],expected)
        self.assertTrue(all(paper.digest(path)==sha for path,sha in before.items()))
        with self.assertRaisesRegex(ValueError,'already exists'):
            paper.analyze(self.args)

    def test_changed_export_and_wrong_full_population_are_rejected(self):
        self.args.allow_subset=False
        with self.assertRaisesRegex(ValueError,'full human'):
            paper.analyze(self.args)
        self.args.allow_subset=True
        path=self.root/'apobec/myoMyo6/rums_mean_by_transcript_id.csv'
        path.write_text(path.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'Changed export'):
            paper.analyze(self.args)
        self.assertFalse(self.args.out_dir.exists())


if __name__=='__main__':
    unittest.main()
