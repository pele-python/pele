import unittest


import networkx as nx
import numpy as np

from pele.landscape import ConnectManager, database2graph
from pele.storage import Database

class TestConnectManager(unittest.TestCase):
    def setUp(self):
        self.db = Database()
        self.nminima = 10
        for i in range(self.nminima):
            e = float(i)
            self.db.addMinimum(e, [e])

    def connect_min(self, m1, m2):
        e = np.random.uniform(0,100)
        self.db.addTransitionState(e, [e], m1, m2)

    def test_gmin(self):
        manager = ConnectManager(self.db, strategy="gmin")
        m0 = self.db.minima()[0]
        
        for i in range(5):
            m1, m2 = manager.get_connect_job()
            self.assertEqual(m1, m0)
            self.connect_min(m1, m2)
    
    def test_random(self):
        manager = ConnectManager(self.db, strategy="random")
        
        for i in range(5):
            m1, m2 = manager.get_connect_job()
            self.connect_min(m1, m2)
        
    def test_combine(self):
        minima = self.db.minima()
        i = 5
        for m1, m2 in zip(minima[:i-1], minima[1:i]):
            self.connect_min(m1, m2)
        for m1, m2 in zip(minima[i:], minima[i+1:]):
            self.connect_min(m1, m2)

        # at this point the minima should be in two disconnected groups
        g = database2graph(self.db)
        self.assertEqual(len(list(nx.connected_components(g))), 2)
            
        manager = ConnectManager(self.db, strategy="combine")
        m1, m2 = manager.get_connect_job()
        self.connect_min(m1, m2)
        
        # they should all be connected now
        g = database2graph(self.db)
        self.assertEqual(len(list(nx.connected_components(g))), 1)
        
    def test_untrap(self):
        # first connect them all randomly
        manager = ConnectManager(self.db, strategy="random")
        while True:
            try:
                m1, m2 = manager.get_connect_job()
            except manager.NoMoreConnectionsError:
                break
            self.connect_min(m1, m2)
        
        for i in range(5):
            try:
                m1, m2 = manager.get_connect_job(strategy="untrap")
            except manager.NoMoreConnectionsError:
                break
            self.connect_min(m1, m2)


    
        
if __name__ == "__main__":
    unittest.main()
        
        
        
        
        
