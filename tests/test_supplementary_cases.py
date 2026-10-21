"""Bounded matched-criterion reconstruction and explicit missing-overlay semantics."""
import json
from pathlib import Path
import sys
import unittest

PAPER=Path(__file__).resolve().parents[1]/'07_manuscript_validation'
sys.path.insert(0,str(PAPER))
import run_supplementary_cases as supplement


class CaseManifestTests(unittest.TestCase):
    def test_reconstruction_is_frozen_and_not_claimed_as_recovered_original(self):
        body=json.loads((PAPER/'supplementary_cases.json').read_text())
        self.assertEqual(body['candidate_transcripts'],985)
        self.assertFalse(body['historical_subset_recovered'])
        self.assertFalse(body['observed_cosmic_variant_overlay_available'])
        self.assertEqual(len({r['transcript_id'] for r in body['s11']}),12)
        self.assertEqual(body['s12']['transcript_name'],'BRCA1-203')


@unittest.skipUnless((PAPER.parent.parent/'sinabro/sinabro/core.py').exists(),'Requires sibling package source')
class MatchedCaseTests(unittest.TestCase):
    def setUp(self):
        supplement.initialize(str(PAPER.parent.parent/'sinabro'))
        self.sequence='C'+'ATG'+'GCTTCA'*4+'TAA'+'G'

    def test_relaxed_paths_share_strict_prefix_and_are_not_shorter(self):
        rows,nodes,edges=supplement.simulate(('s11','toy',self.sequence,0,10,42,1000))
        self.assertEqual(len(rows),10)
        self.assertTrue(all(r['blosum_neutral_steps']>=r['identity_neutral_steps'] for r in rows))
        self.assertEqual(nodes,[])
        self.assertEqual(edges,[])

    def test_simulation_graph_has_real_edges_and_no_invented_cosmic_annotations(self):
        rows,nodes,edges=supplement.simulate(('s12','toy',self.sequence,0,3,42,1000))
        self.assertEqual(len(rows),3)
        ids={r['node'] for r in nodes}
        self.assertTrue(all(e['source'] in ids and e['target'] in ids for e in edges))
        self.assertTrue(all(r['observed_cosmic_status']=='not_available' for r in nodes))
        self.assertTrue(any(r['hamming_distance']==0 for r in nodes))

    def test_incomplete_blosum_path_is_not_reported_as_success(self):
        with self.assertRaisesRegex(ValueError,'Incomplete BLOSUM path'):
            supplement.simulate(('s11','toy',self.sequence,0,1,42,1))


if __name__=='__main__':
    unittest.main()
