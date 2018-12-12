import unittest
import numpy as np

from pele.storage import Database
from pele.landscape import ConnectManager
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph, TreeLeastCommonAncestor, Tree

_show = False

def create_random_database(nmin=20, nts=20):
    db = Database()
    
    for i in range(nmin):
        energy = np.random.uniform(-1., 10.)
        x = [energy]
        db.addMinimum(energy, x)
    
    manager = ConnectManager(db, verbosity=0)
    for i in range(nts//2):
        m1, m2 = manager.get_connect_job("gmin")
        energy = max([m1.energy, m2.energy]) + np.random.uniform(1,10)
        x = [energy]
        db.addTransitionState(energy, x, m1, m2)
    for i in range(nts//2):
        m1, m2 = manager.get_connect_job("random")
        energy = max([m1.energy, m2.energy]) + np.random.uniform(1,5)
        x = [energy]
        db.addTransitionState(energy, x, m1, m2)
    
    return db
    


class TestDisconnectivityGraph(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.db = create_random_database(20, 60)
        self.tsgraph = database2graph(self.db)

    def test_basic(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        dgraph.plot()
        layout = dgraph.get_tree_layout()
        
        dgraph.draw_minima(self.db.minima()[:2])

    def test_show_minima(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        dgraph.plot(show_minima=True)

    def test_color_groups(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        groups = [ self.db.minima()[:5],
                  self.db.minima()[5:10],
                  self.db.minima()[10:15],
                  self.db.minima()[15:18],
                  ]
        dgraph.color_by_group(groups)
        dgraph.plot(linewidth=2.)
        if _show:
            from matplotlib import pyplot as plt
            plt.title("color by group")
#             dgraph.show()
    
    def test_color_groups_many(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        groups = []
        for m in self.db.minima():
            groups.append([m])
            if len(groups) > 13:
                break
        dgraph.color_by_group(groups)
        dgraph.plot(linewidth=2.)
        if _show:
            from matplotlib import pyplot as plt
            plt.title("color by group")
#             dgraph.show()

    def test_color_value(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        def minimum_to_value(m):
            if m.energy < 5.:
                return m.energy
            else:
                return None
        dgraph.color_by_value(minimum_to_value)
        dgraph.plot(linewidth=2.)
        dgraph.label_minima({self.db.minima()[0]: "gmin", 
                             self.db.minima()[1]: "2nd lowest"})
        if _show:
            from matplotlib import pyplot as plt
            plt.title("color by value")
            dgraph.show()

    def test_Emax(self):
        emax = self.db._highest_energy_minimum().energy
        m1 = self.db.addMinimum(emax-5.1, [0.])
        m2 = self.db.addMinimum(emax-5.2, [0.])
        self.db.addTransitionState(emax-5, [0.], m1, m2)
        self.db.addMinimum(emax-5.3, [0.]) # so this will need to be removed
        dgraph = DisconnectivityGraph(self.tsgraph, Emax=emax-2, subgraph_size=1, order_by_energy=True)
        dgraph.calculate()
        dgraph.plot()
    
    def test_gmin_not_connected(self):
        m0 = self.db.get_lowest_energy_minimum()
        self.db.addMinimum(m0.energy - 10., m0.coords)
        
        tsgraph = database2graph(self.db)
        dgraph = DisconnectivityGraph(tsgraph, include_gmin=True)
        dgraph.calculate()
        dgraph.plot()
    
    def test_2_minima(self):
        db = create_random_database(nmin=3, nts=1)
        
        tsgraph = database2graph(self.db)
        dgraph = DisconnectivityGraph(tsgraph, include_gmin=True)
        dgraph.calculate()
        dgraph.plot()


class TestTreeLeastCommonAncestor(unittest.TestCase):
    def test(self):
        t1 = Tree()
        t2_1 = Tree(t1)
        t2_2 = Tree(t1)
        lca = TreeLeastCommonAncestor([t2_1, t2_2])
        self.assertEqual(lca.least_common_ancestor, t1)
        
        t3_1 = Tree(t2_1)
        t3_2 = Tree(t2_1)
        t3_3 = Tree(t2_2)
        
        self.assertEqual(TreeLeastCommonAncestor([t3_1, t3_2]).least_common_ancestor, t2_1)
        self.assertEqual(TreeLeastCommonAncestor([t3_1, t3_2, t3_3]).least_common_ancestor, t1)
        self.assertEqual(TreeLeastCommonAncestor([t2_1, t3_1]).least_common_ancestor, t2_1)
        self.assertEqual(TreeLeastCommonAncestor([t2_2, t3_1]).least_common_ancestor, t1)
        
        lca = TreeLeastCommonAncestor([t3_1, t3_2])

        try:
            self.assertItemsEqual(lca.get_all_paths_to_common_ancestor(), [t2_1, t3_1, t3_2])
        except AttributeError:
            self.assertCountEqual(lca.get_all_paths_to_common_ancestor(), [t2_1, t3_1, t3_2])

    def test_make_branch(self):
        t1 = Tree()
        t2 = t1.make_branch()
        self.assertEqual(t2.parent, t1)
        self.assertEqual(t1.number_of_subtrees(), 2)

if __name__ == "__main__":
    unittest.main()

