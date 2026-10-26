"""Bounded matched-criterion reconstruction and explicit missing-overlay semantics."""
import json
import copy
from pathlib import Path
import sys
import unittest

PAPER=Path(__file__).resolve().parents[1]/'07_robustness_analysis'
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

    def test_graph_audit_reconstructs_dna_and_checks_every_complete_path(self):
        rows,nodes,edges=supplement.simulate(('s12','toy',self.sequence,0,8,42,1000))
        audit=supplement.audit_graph(rows,{r['node']:r for r in nodes},edges,self.sequence)
        self.assertTrue(audit['graph_audit_passed'])
        self.assertEqual(audit['complete_paths'],8)
        self.assertEqual(audit['single_substitution_edges_checked'],
                         sum(r['blosum_neutral_steps']+1 for r in rows))

    def test_graph_audit_rejects_changed_coordinate_and_chronology(self):
        rows,nodes,edges=supplement.simulate(('s12','toy',self.sequence,0,3,42,1000))
        altered=copy.deepcopy({r['node']:r for r in nodes})
        next(iter(altered.values()))['cds_substitutions']='1C>A'
        with self.assertRaisesRegex(ValueError,'coordinate/reference'):
            supplement.audit_graph(rows,altered,edges,self.sequence)
        altered_edges=copy.deepcopy(edges);altered_edges[0]['step']=2
        with self.assertRaisesRegex(ValueError,'Missing or duplicate graph step'):
            supplement.audit_graph(rows,{r['node']:r for r in nodes},altered_edges,self.sequence)

    def test_path_seeds_and_graph_states_do_not_depend_on_work_blocks(self):
        whole=supplement.simulate(('s12','toy',self.sequence,0,7,42,1000))
        parts=[supplement.simulate(('s12','toy',self.sequence,start,count,42,1000))
               for start,count in [(0,3),(3,4)]]
        self.assertEqual(whole[0],parts[0][0]+parts[1][0])
        self.assertEqual({r['node']:r for r in whole[1]},
                         {r['node']:r for part in parts for r in part[1]})
        self.assertEqual(whole[2],parts[0][2]+parts[1][2])


if __name__=='__main__':
    unittest.main()
