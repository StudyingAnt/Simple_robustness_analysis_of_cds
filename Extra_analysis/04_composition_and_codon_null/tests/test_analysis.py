"""Numerical-identifiability, held-out data isolation and stochastic checks."""
import hashlib
import os
from pathlib import Path
import sys
import unittest

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import analyze_composition as composition
import run_codon_null as null
from common import np


class RegressionTests(unittest.TestCase):
    def test_redundant_covariates_recover_known_response(self):
        rng=np.random.default_rng(105)
        x=rng.normal(size=(80,3))
        y=3+2*x[:,0]-.5*x[:,1]+.7*x[:,2]
        redundant=np.column_stack([x,x[:,0]+x[:,1],x[:,2],np.ones(80)])
        prediction,rank=composition.fit_predict(redundant,y,np.arange(60),np.arange(60,80))
        np.testing.assert_allclose(prediction,y[60:],atol=1e-10)
        self.assertEqual(rank,4)

    def test_test_responses_do_not_enter_training(self):
        rng=np.random.default_rng(106)
        x=rng.normal(size=(80,5));y=rng.normal(size=80)
        a,_=composition.fit_predict(x,y,np.arange(60),np.arange(60,80))
        y[60:]+=10000
        b,_=composition.fit_predict(x,y,np.arange(60),np.arange(60,80))
        np.testing.assert_array_equal(a,b)


class NullTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=Path(os.environ.get('SINABRO_SOURCE',str(HERE.parents[2]/'sinabro')))
        package=source/'sinabro'
        h=hashlib.sha256()
        for path in sorted(package.rglob('*.py')):
            h.update(str(path.relative_to(source)).encode()+b'\0');h.update(path.read_bytes())
        null.initialize(str(source),h.hexdigest())

    def test_recoding_protein_and_first_step_probability(self):
        sequence='A'+'ATGGCTTTTTGA'+'C'
        success=np.zeros(5);expect=np.zeros(5);variance=np.zeros(5)
        hashes=set()
        for variant in range(11):
            counts,meta=null.compute(('known_protein',sequence,variant,40,105,10000))
            self.assertTrue(meta['protein_preserved'])
            self.assertEqual(meta['protein_sha256'],hashlib.sha256(b'MAF*').hexdigest())
            self.assertEqual(meta['cds_nt'],12)
            if variant==0:self.assertEqual(meta['sequence_sha256'],hashlib.sha256(sequence.encode()).hexdigest())
            hashes.add(meta['sequence_sha256'])
            p=np.array(meta['weighted_one_step'])
            success+=(counts>=1).sum(axis=0);expect+=40*p;variance+=40*p*(1-p)
        self.assertGreater(len(hashes),4)
        self.assertTrue((np.abs((success-expect)/np.sqrt(variance))<6).all())

    def test_exact_replay_and_nested_path_prefix(self):
        job=('replay','C'+'ATGGCTGAATTCTAA'+'G',2,12,20261025,10000)
        a,ma=null.compute(job);b,mb=null.compute(job)
        np.testing.assert_array_equal(a,b)
        self.assertEqual(ma['sequence_sha256'],mb['sequence_sha256'])
        larger=list(job);larger[3]=18
        c,_=null.compute(tuple(larger))
        np.testing.assert_array_equal(a,c[:12])


if __name__=='__main__':unittest.main()
