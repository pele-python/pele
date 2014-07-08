import unittest
import numpy as np

from pele.storage import Database
from pele.landscape import ConnectManager
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph

_show = False

def create_random_database(nmin=20, nts=20):
    db = Database()
    
    for i in xrange(nmin):
        energy = np.random.uniform(-1., 10.)
        x = [energy]
        db.addMinimum(energy, x)
    
    manager = ConnectManager(db, verbosity=0)
    for i in xrange(nts/2):
        m1, m2 = manager.get_connect_job("gmin")
        energy = max([m1.energy, m2.energy]) + np.random.uniform(1,10)
        x = [energy]
        db.addTransitionState(energy, x, m1, m2)
    for i in xrange(nts/2):
        m1, m2 = manager.get_connect_job("random")
        energy = max([m1.energy, m2.energy]) + np.random.uniform(1,5)
        x = [energy]
        db.addTransitionState(energy, x, m1, m2)
    
    return db
    


class TestDisconnectivityGraph(unittest.TestCase):
    def setUp(self):
#         np.random.seed(1)
        self.db = create_random_database(20, 60)
        self.tsgraph = database2graph(self.db)

    def test_basic(self):
        dgraph = DisconnectivityGraph(self.tsgraph)
        dgraph.calculate()
        dgraph.plot()

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


if __name__ == "__main__":
    unittest.main()
