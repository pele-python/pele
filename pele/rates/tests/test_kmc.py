import unittest
import numpy as np

from pele.rates._kmc import KineticMonteCarlo
from pele.rates._rate_calculations import GraphReduction
from test_graph_transformation import _MakeRandomGraph, _three_state_graph


class TestKMC(unittest.TestCase):
    def setUp(self):
        graph = _three_state_graph()
        self.kmc = KineticMonteCarlo(graph)
    
    def test_mfp(self):
        time = self.kmc.mean_first_passage_time(0, [1], niter=1000)
        self.assertAlmostEqual(time, 1.0, delta=.1)
    

class TestKMC_GraphReduction(unittest.TestCase):
    def compare(self, A, B, nnodes=10, nedges=20, weights=None, x=1):
        print ""
        maker = _MakeRandomGraph(nnodes=nnodes, nedges=nedges, node_set=A+B+[x])
        graph = maker.run()
        graph_backup = graph.copy()
        reducer = GraphReduction(graph, A, B, weights=weights)
        kmc = KineticMonteCarlo(graph_backup, debug=False)
        
        # test compute_committor_probability()
        PxB = reducer.compute_committor_probability(x)
        PxB_kmc = kmc.committor_probability(x, A, B, niter=1000)
        print "committor probability    ", x, "->", B, "=", PxB
        print "committor probability kmc", x, "->", B, "=", PxB_kmc
        self.assertAlmostEqual(PxB, PxB_kmc, delta=.1)
        
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
         
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
        
        paB = kmc.committor_probability(A[0], [A[0]], B, niter=1000)
        print "the committor probability a->B", paB
        print "graph reduction committor prob", reducer.get_committor_probabilityAB(A[0])
        self.assertAlmostEqual(paB, reducer.get_committor_probabilityAB(A[0]), delta=.1)

    def test(self):
        A = [0]
        B = [1]
        self.compare(A, B, x=3)

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
        self.compare(A, B, weights=weights, x=4)


    def test_big_group(self, nnodes=10, nedges=20):
        A = range(8)
        B = [9]
        self.compare(A, B)
    
    def test_committor_probabilities(self, nnodes=10, nedges=20):
        A = [0,1,2,3]
        B = [8,9]
        xx = 5
        maker = _MakeRandomGraph(nnodes=nnodes, nedges=nedges, node_set=A+B)
        graph = maker.run()
        graph_backup = graph.copy()
        kmc = KineticMonteCarlo(graph_backup, debug=False)
        reducer = GraphReduction(graph, A, B)
        
        nodes = set(A + B + [xx])
        PxB = reducer.compute_committor_probabilities(nodes)
        for x in nodes:
            self.assertIn(x, PxB)
            
        for x in nodes:
            PxB_kmc = kmc.committor_probability(x, A, B, niter=1000)
            self.assertAlmostEqual(PxB[x], PxB_kmc, delta=.1)
  


if __name__ == "__main__":
    unittest.main()
