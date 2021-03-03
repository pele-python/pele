import unittest

import networkx as nx

from pele.landscape import TSGraph


def create_random_database(nmin=20, nts=None, natoms=2):
    """
    create a database for test purposes
    """
    from pele.storage import Database
    import numpy as np

    if nts is None:
        nts = nmin
    db = Database()
    # generate random structures
    minlist = []
    for i in range(nmin):
        coords = np.random.uniform(-1, 1, natoms * 3)
        e = float(i)  # make up a fake energy
        minlist.append(db.addMinimum(e, coords))
    # add random transition states
    for i in range(nts):
        j1, j2 = 1, 1
        while j1 == j2:
            j1, j2 = np.random.randint(0, nmin, 2)
        m1, m2 = minlist[j1], minlist[j2]
        coords = np.random.uniform(-1, 1, natoms * 3)
        e = float(j1 + j2)
        db.addTransitionState(e, coords, m1, m2)
    return db


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.db = create_random_database()

    def test_hi(self):
        graph = TSGraph(self.db)
        ng = graph.graph.number_of_nodes()
        nd = len(self.db.minima())
        self.assertEqual(ng, nd, "all nodes not imported")

        ntsg = graph.graph.number_of_edges()
        ntsd = len(self.db.transition_states())
        self.assertEqual(ntsg, ntsd, "all transition states not imported")

    def test_minima_keyword(self):
        nmin = 3
        minima = list(self.db.minima())
        minima = minima[:nmin]
        graph = TSGraph(self.db, minima=minima, no_edges=True)
        ng = graph.graph.number_of_nodes()
        self.assertEqual(ng, nmin)

        self.assertEqual(graph.graph.number_of_edges(), 0)

    def test_merge_minima(self):
        graph = TSGraph(self.db)
        # select two minima with at least one edge
        minima = []
        for n in graph.graph.nodes():
            if graph.graph.degree(n) > 0:
                minima.append(n)
            if len(minima) == 2:
                break
        nedges = [graph.graph.degree(n) for n in minima]
        min1, min2 = minima
        neighbors = graph.graph.neighbors(min2)
        graph.mergeMinima(min1, min2)
        self.db.mergeMinima(min1, min2)

        self.assertNotIn(min2, graph.graph.nodes())

        # make sure the edges were copied
        for n in neighbors:
            if n == min1: continue
            edge_exists = (graph.graph.has_edge(min1, n)
                           or graph.graph.has_edge(n, min1))
            self.assert_(edge_exists)


        # make sure the new edges have transition state data attached
        for v in graph.graph.neighbors(min1):
            data = graph.graph.get_edge_data(v, min1)
            self.assertIn("ts", data.keys())

        # make sure min2 is not in database
        self.assertNotIn(min2, self.db.minima())

    def test_networkx(self):
        """check how networkx works"""
        graph = nx.Graph()
        graph.add_edge(1, 2, ts=1)
        graph.add_edge(2, 1, ts=2)
        self.assertEqual(graph.number_of_edges(), 1)

        data = graph.get_edge_data(1, 2)
        self.assertIsNotNone(data)
        self.assertDictEqual(graph.get_edge_data(1, 2),
                             graph.get_edge_data(2, 1))

        self.assertTrue(graph.has_edge(1, 2))
        self.assertTrue(graph.has_edge(2, 1))

    def test_connected_components(self):
        tsgraph = TSGraph(self.db)
        cc = nx.connected_components(tsgraph.graph)
        # networkx changed the function so now cc is an iterator over sets
        cc = [list(c) for c in cc] 
        for nodes in cc:
            for u, v in zip(nodes[:-1], nodes[1:]):
                self.assertTrue(tsgraph.areConnected(u, v))

        for nodes1, nodes2 in zip(cc[:-1], cc[1:]):
            for u in nodes1:
                for v in nodes2:
                    self.assertFalse(tsgraph.areConnected(u, v))


if __name__ == "__main__":
    unittest.main()

