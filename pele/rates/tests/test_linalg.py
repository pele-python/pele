from __future__ import print_function
from __future__ import absolute_import
import unittest
from .test_graph_transformation import _three_state_rates, _MakeRandomGraph

from pele.rates._rates_linalg import CommittorLinalg, MfptLinalgSparse, TwoStateRates

class TestLinalg3(unittest.TestCase):
    def setUp(self):
        self.rates = _three_state_rates()
        # all rates after graph renormalization should be 1.0
        self.final_rate = 1.0

    def _test_rate(self, i, j):
        reducer = TwoStateRates(self.rates, [i], [j])
        reducer.compute_rates()
        rAB = reducer.get_rate_AB() 
        self.assertAlmostEqual(rAB, self.final_rate, 7)
        
        reducer.compute_committors()
        rAB_ss = reducer.get_rate_AB_SS()
        print("kSS", rAB_ss)
        self.assertAlmostEqual(rAB_ss, 1.5, 7)
        
    def test01(self):
        self._test_rate(0,1)

    def test12(self):
        self._test_rate(1,2)

    def test02(self):
        self._test_rate(0,2)

class TestNgtCpp10(unittest.TestCase):
    def setUp(self):
        from .test_cpp_ngt import make_rates_complete
        self.rates = make_rates_complete(nnodes=10)
        # all rates after graph renormalization should be 1.0
        self.true_kAB = 5.1013138820442565
        self.true_kBA = 3
        self.true_kAB_SS = 19.933950145409426
        self.true_kBA_SS = 6.970856553435547

    def _test_rate(self, A, B):
        reducer = TwoStateRates(self.rates, A, B)
        reducer.compute_rates()
        reducer.compute_committors()
        kAB = reducer.get_rate_AB()
        self.assertAlmostEqual(kAB, self.true_kAB, 7)
#        kBA = reducer.get_rate_BA()
#        self.assertAlmostEqual(kBA, self.true_kBA, 7)

        rAB_SS = reducer.get_rate_AB_SS()
        self.assertAlmostEqual(rAB_SS, self.true_kAB_SS, 7)
#        rBA_SS = reducer.get_rate_BA_SS()
#        self.assertAlmostEqual(rBA_SS, self.true_kBA_SS, 7)

    def test01(self):
        self._test_rate([0,1,2],[3,4,5])


class TestLinalgRandom(unittest.TestCase):
    def do_check(self, A, B, nnodes=20, nedges=20):
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        rates = maker.make_rates()
        reducer = TwoStateRates(rates, A, B)
        reducer.compute_rates()
        reducer.compute_committors()
        
        from pele.rates._ngt_cpp import NGT
        ngt = NGT(rates, A, B)
        ngt.compute_rates()
        
        self.assertAlmostEqual(reducer.get_rate_AB(), ngt.get_rate_AB(), 7)
        self.assertAlmostEqual(reducer.get_rate_AB_SS(), ngt.get_rate_AB_SS(), 7)


    def test(self):
        A, B = [0], [1]
        self.do_check(A, B)
 
    def test_setA(self):
        A, B = [0, 1, 2], [3]
        self.do_check(A, B)
  
    def test_setAB(self):
        A, B = [0, 1, 2], [3, 4, 5, 6]
        self.do_check(A, B)

    def test_weakly_connected(self):
        rates = {}
        nnodes = 5
        for n in range(nnodes-1):
            rates[(n,n+1)] = 1.#np.random.rand()           
            rates[(n+1,n)] = 1.#np.random.rand()
        # add a disconnected edge
#         rates[(nnodes, nnodes+1)] = 1.
#         rates[(nnodes+1, nnodes)] = 1.
        
        for a in range(nnodes-1):
            calc = MfptLinalgSparse(rates, [a])
            times = calc.compute_mfpt()
            self.assertGreater(min(times.values()), 0)


if __name__ == "__main__":
    unittest.main()

