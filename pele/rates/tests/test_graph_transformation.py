import unittest
import numpy as np
import networkx as nx

from pele.rates import RateCalculation
from pele.rates._rate_calculations import GraphReduction, kmcgraph_from_rates

np.random.seed(0)

def make_rates_complete(nnodes=10):
    rates = dict()
    for i in range(nnodes):
        for j in range(i+1,nnodes):
            rates[(i,j)] = float(i+j) / (i+1)
            rates[(j,i)] = float(i+j) / (j+1)
    return rates


class _MakeRandomGraph(object):
    def __init__(self, nnodes=10, nedges=20, node_set=None):
        self.nodes = np.array(list(range(nnodes)))
        self.nedges = nedges
        self.node_set = set(node_set)
        self.rates = dict()
        self.graph = nx.Graph()
        self.graph.add_nodes_from(self.nodes)
    
    def node_set_connected(self):
        u = next(iter(self.node_set))
        cc = nx.node_connected_component(self.graph, u)
        cc = set(cc)
        return len(cc.intersection(self.node_set)) == len(self.node_set)
    
    def add_random_edge(self):
        np.random.shuffle(self.nodes)
        u, v = self.nodes[:2]
        r1, r2 = np.random.uniform(.5,1.5, 2)
        self.rates[(u,v)] = r1
        self.rates[(v,u)] = r2
        self.graph.add_edge(u,v)
#         print len(self.rates) / 2, "edges"
        return u, v

    def make_rates(self):
        nnodes = len(self.nodes)
        if self.node_set is not None:
            # add edges until u and v are connected    
            while not self.node_set_connected():
                self.add_random_edge()
        nedges = min(self.nedges, nnodes*(nnodes-1))
        while len(self.rates) < nedges:
            self.add_random_edge()
#        print "made random graph with", len(self.nodes), "nodes and", len(self.rates) / 2, "edges"
        return self.rates
    
    def run(self):
        self.make_rates()
        return kmcgraph_from_rates(self.rates)
    
    


def _three_state_rates():
    tmatrix = [ [0., 1., 1.,], [1., 0., 1.,], [1., 1., 0.] ]
    rates = dict()
    for i in range(3):
        for j in range(3):
            if i != j:
                rates[(i,j)] = tmatrix[i][j]
    return rates

def _three_state_graph():
    return kmcgraph_from_rates(_three_state_rates())

class TestGraphReduction3(unittest.TestCase):
    def setUp(self):
        self.rates = _three_state_rates()
        # all rates after graph renormalization should be 1.0
        self.final_rate = 1.0

    def _test_rate(self, i, j):
        reducer = GraphReduction(self.rates, [i], [j], debug=True)
        reducer.check_graph()
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
        reducer.check_graph()
        self.assertEqual(reducer.graph.number_of_nodes(), 2)
        self.assertEqual(reducer.graph.number_of_edges(), 4)
        self.assertAlmostEqual(rAB, self.final_rate, 7)
        self.assertAlmostEqual(rBA, self.final_rate, 7)

    def test01(self):
        self._test_rate(0,1)

    def test12(self):
        self._test_rate(1,2)

    def test02(self):
        self._test_rate(0,2)

class TestGraphReductionRandom(unittest.TestCase):
    def do_check(self, A, B, nnodes=20, nedges=20):
        maker = _MakeRandomGraph(nnodes=20, nedges=20, node_set=A+B)
        graph = maker.run()
        reducer = GraphReduction(maker.rates, A, B, debug=False)  
        reducer.check_graph()
        reducer.compute_rates()
        rAB = reducer.get_rate_AB()
        rBA = reducer.get_rate_BA()
        reducer.check_graph()
        self.assertEqual(reducer.graph.number_of_nodes(), len(A) + len(B))
        if len(A) == 1 and len(B) == 1:
            self.assertLessEqual(reducer.graph.number_of_edges(), 4)
             
    def test(self):
        A, B = [0], [1]
        self.do_check(A, B)
 
    def test_setA(self):
        A, B = [0, 1, 2], [3]
        self.do_check(A, B)
  
    def test_setAB(self):
        A, B = [0, 1, 2], [3, 4, 5, 6]
        self.do_check(A, B)


if __name__ == "__main__":
    unittest.main()

