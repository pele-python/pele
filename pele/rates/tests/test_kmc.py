import unittest
import numpy as np

from pele.rates._kmc import KineticMonteCarlo
from pele.rates._rate_calculations import GraphReduction
from test_graph_transformation import _MakeRandomGraph


class TestKMC(unittest.TestCase):
    def setUp(self):
        from pele.rates.tests.test_graph_transformation import _three_state_graph
        graph = _three_state_graph()
        self.kmc = KineticMonteCarlo(graph)
    
    def test_mfp(self):
        time = self.kmc.mean_first_passage_time(0, [1], niter=1000)
        self.assertAlmostEqual(time, 1.0, delta=.1)
    

class TestKMC_GraphReduction(unittest.TestCase):
    def compare(self, A, B, nnodes=10, nedges=20, weights=None):
        maker = _MakeRandomGraph(nnodes=nnodes, nedges=nedges, node_set=A+B)
        graph = maker.run()
        graph_backup = graph.copy()
        reducer = GraphReduction(graph, A, B, weights=weights)  
        rAB, rBA = reducer.compute_rates()
         
        kmc = KineticMonteCarlo(graph_backup, debug=False)
        rAB_KMC = kmc.mean_rate(A, B, niter=1000, weights=weights)
         
        print "NGT rate A->B", rAB
        print "KMC rate A->B", rAB_KMC
        print "normalized difference", (rAB - rAB_KMC)/rAB 
        self.assertLess(abs(rAB - rAB_KMC)/rAB, .1)
         
        rBA_KMC = kmc.mean_rate(B, A, niter=1000, weights=weights)
         
        print "NGT rate B->A", rBA
        print "KMC rate B->A", rBA_KMC
        print "normalized difference", (rBA - rBA_KMC)/rBA
        self.assertLess(abs(rBA - rBA_KMC)/rBA, .1)
    
    def test(self):
        A = [0]
        B = [1]
        self.compare(A, B)

    def test_group(self, nnodes=10, nedges=20):
        A = [0,1]
        B = [2, 3]
        self.compare(A, B)

    def test_weights(self, nnodes=10, nedges=20):
        A = [0,1]
        B = [2, 3]
        weights = dict()
        for x in A+B:
            weights[x] = np.random.uniform(0,1)
        self.compare(A, B, weights=weights)


if __name__ == "__main__":
    unittest.main()
