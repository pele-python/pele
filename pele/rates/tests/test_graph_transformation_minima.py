import unittest
import os
import numpy as np
import networkx as nx

from pele.systems import LJCluster
from pele.landscape import ConnectManager
from pele.thermodynamics import get_thermodynamic_information
from pele.utils.disconnectivity_graph import database2graph

from pele.rates import RateCalculation
from pele.rates._rate_calculations import GraphReduction, kmcgraph_from_rates

from test_graph_transformation import _MakeRandomGraph

def create_random_database(system, db, nmin=20, nts=10):
    if db.number_of_minima() < nmin:
        bh = system.get_basinhopping(db, outstream=None)
        bh.run(50)
    
    if db.number_of_transition_states() < nts:
        manager = ConnectManager(db, strategy="gmin")
        for i in range(nts):
            min1, min2 = manager.get_connect_job("gmin")
            connect = system.get_double_ended_connect(min1, min2, db, verbosity=0)
            connect.connect()
        

    connect = system.get_double_ended_connect(db.minima()[0], db.minima()[-1], db, verbosity=0)
    connect.connect()
        
    get_thermodynamic_information(system, db, nproc=2)
    return db

class TestGraphRatesLJ(unittest.TestCase):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        dbfname = current_dir + "/lj13.db"
        self.system = LJCluster(13)
        self.system.params.structural_quench_params.tol = 1e-6
        self.db = self.system.create_database(dbfname)
        self.ndof = self.system.get_ndof()
        
        create_random_database(self.system, self.db, 10, 20)
    
    def test(self):
        A = [self.db.minima()[0]]
        B = [self.db.minima()[-1]]
        
        graph = database2graph(self.db)
        rcalc = RateCalculation(graph, A, B, ndof=self.ndof, T=1.)
        rAB, rBA = rcalc.compute_rates()
        print "rates", rAB, rBA

    def test2(self):
        A = self.db.minima()[:2]
        B = self.db.minima()[2:4]

        graph = database2graph(self.db)
        rcalc = RateCalculation(graph, A, B, T=1., ndof=self.ndof)
        rAB, rBA = rcalc.compute_rates()
        print "rates", rAB, rBA

if __name__ == "__main__":
    unittest.main()