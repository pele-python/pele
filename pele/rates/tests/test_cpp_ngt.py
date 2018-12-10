from __future__ import absolute_import
import unittest
import numpy as np
import networkx as nx

from pele.rates._ngt_cpp import NGT
from .test_graph_transformation import _MakeRandomGraph, _three_state_rates, make_rates_complete



class TestNgtCpp3(unittest.TestCase):
    def setUp(self):
        self.rates = _three_state_rates()
        # all rates after graph renormalization should be 1.0
        self.final_rate = 1.0
        self.final_rate_SS = 1.5

    def _test_rate(self, i, j):
        reducer = NGT(self.rates, [i], [j], debug=False)
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
        self.assertAlmostEqual(rAB, self.final_rate, 7)
        self.assertAlmostEqual(rBA, self.final_rate, 7)

        rAB_SS = reducer.get_rate_AB_SS()
        rBA_SS = reducer.get_rate_BA_SS()
        self.assertAlmostEqual(rAB_SS, self.final_rate_SS, 7)
        self.assertAlmostEqual(rBA_SS, self.final_rate_SS, 7)

    def test01(self):
        self._test_rate(0,1)

    def test12(self):
        self._test_rate(1,2)

    def test02(self):
        self._test_rate(0,2)

class TestNgtCpp10(unittest.TestCase):
    def setUp(self):
        self.rates = make_rates_complete(nnodes=10)
        # all rates after graph renormalization should be 1.0
        self.true_kAB = 5.1013138820442565
        self.true_kBA = 3
        self.true_kAB_SS = 19.933950145409426
        self.true_kBA_SS = 6.970856553435547

    def _test_rate(self, A, B):
        reducer = NGT(self.rates, A, B, debug=False)
        reducer.compute_rates()
        kAB = reducer.get_rate_AB()
        kBA = reducer.get_rate_BA()
        self.assertAlmostEqual(kAB, self.true_kAB, 7)
        self.assertAlmostEqual(kBA, self.true_kBA, 7)

        rAB_SS = reducer.get_rate_AB_SS()
        rBA_SS = reducer.get_rate_BA_SS()
        self.assertAlmostEqual(rAB_SS, self.true_kAB_SS, 7)
        self.assertAlmostEqual(rBA_SS, self.true_kBA_SS, 7)

    def test01(self):
        self._test_rate([0,1,2],[3,4,5])


class TestNgtCppRandom(unittest.TestCase):
    def do_check(self, A, B, nnodes=20, nedges=20):
        np.random.seed(0)
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        maker.run()
        reducer = NGT(maker.rates, A, B, debug=False)  
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
        
    def compare_linalg(self, A, B, nnodes=20, nedges=20):
        from pele.rates._rates_linalg import TwoStateRates
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        maker.run()
        
        reducer = NGT(maker.rates, A, B)
        reducer.compute_rates_and_committors()
        committors = reducer.get_committors()
        
        la = TwoStateRates(maker.rates, A, B)
#        la.compute_rates()
        la.compute_committors()
        qla = la.committor_dict
        for n, qla in la.committor_dict.items():
            self.assertAlmostEqual(qla, committors[n], 7)
        
    
    def test(self):
        A, B = [0], [1]
        self.do_check(A, B)
        self.compare_linalg(A, B)
 
    def test_setA(self):
        A, B = [0, 1, 2], [3]
        self.do_check(A, B)
        self.compare_linalg(A, B)
  
    def test_setAB(self):
        A, B = [0, 1, 2], [3, 4, 5, 6]
        self.do_check(A, B)
        self.compare_linalg(A, B)


if __name__ == "__main__":
    unittest.main()

