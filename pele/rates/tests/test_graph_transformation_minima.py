from __future__ import print_function
from __future__ import absolute_import
import unittest
import os
import numpy as np
import networkx as nx

from pele.systems import LJCluster
from pele.landscape import ConnectManager
from pele.thermodynamics import get_thermodynamic_information
from pele.utils.disconnectivity_graph import database2graph

from pele.rates import RateCalculation, RatesLinalg
from pele.rates._rate_calculations import GraphReduction, kmcgraph_from_rates

from .test_graph_transformation import _MakeRandomGraph

def create_random_database(system, db, nmin=20, nts=10):
    if db.number_of_minima() < nmin:
        bh = system.get_basinhopping(db, outstream=None)
        bh.run(50)
    
    if db.number_of_transition_states() < nts:
        manager = ConnectManager(db, strategy="gmin")
        for i in range(nts):
            try:
                min1, min2 = manager.get_connect_job("gmin")
            except Exception as e:
                if not "couldn't find any random minima pair to connect" in str(e):
                    raise
                
                    
            connect = system.get_double_ended_connect(min1, min2, db, verbosity=0)
            connect.connect()
        

    connect = system.get_double_ended_connect(db.minima()[0], db.minima()[-1], db, verbosity=0)
    connect.connect()
        
    get_thermodynamic_information(system, db, nproc=2)
    return db

class TestGraphRatesLJ(unittest.TestCase):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        dbfname = os.path.join(current_dir, "lj15.sqlite")
        print(dbfname)
        self.system = LJCluster(15)
        self.system.params.structural_quench_params.tol = 1e-6
        self.db = self.system.create_database(dbfname, createdb=False)
        
#        create_random_database(self.system, self.db, 10, 20)
#        get_thermodynamic_information(self.system, self.db, nproc=2)
    
    def do_check_rates(self, A, B):
        rcalc = RateCalculation(self.db.transition_states(), A, B)
        rcalc.compute_rates()
        rAB = rcalc.get_rate_AB()
        rBA = rcalc.get_rate_BA()
        
        rla = RatesLinalg(self.db.transition_states(), A, B)
        rAB_la = rla.compute_rates()
        
        self.assertAlmostEqual(rAB, rAB_la, 7)
        
    def do_check_committors(self, A, B):
        rcalc = RateCalculation(self.db.transition_states(), A, B)
        rcalc.compute_rates_and_committors()
        committors = rcalc.get_committors()
        
        rla = RatesLinalg(self.db.transition_states(), A, B)
        cla = rla.compute_committors()
        
        for m, qla in cla.items():
            self.assertAlmostEqual(committors[m], qla, 7)
        
    def test(self):
        A = [self.db.minima()[0]]
        B = [self.db.minima()[-1]]
        self.do_check_rates(A, B)
        self.do_check_committors(A, B)

    def test2(self):
        A = self.db.minima()[:2]
        B = self.db.minima()[2:4]
        self.do_check_rates(A, B)
        self.do_check_committors(A, B)

class TestOptimCollagen(unittest.TestCase):
    """test a known value for a large database"""
    def setUp(self):
        from pele.utils.optim_compatibility import OptimDBConverter
        from pele.storage import Database
        ndof = 10 # wrong, but who cares.
        self.db = Database()
        current_dir = os.path.dirname(__file__)
        converter = OptimDBConverter(self.db, ndof=ndof, mindata=current_dir+"/collagen.min.data", 
                                     tsdata=current_dir+"/collagen.ts.data", assert_coords=False)
        converter.convert_no_coords()
    
    def test1(self):
        m1 = self.db.getMinimum(1)
        m2 = self.db.getMinimum(2)
        m3 = self.db.getMinimum(3)
        m4 = self.db.getMinimum(4)
        
        rcalc = RateCalculation(self.db.transition_states(), [m1], [m2], T=0.592)
        rcalc.compute_rates()
        self.assertAlmostEqual(rcalc.get_rate_AB(), 7106337458., delta=1e4)
        self.assertAlmostEqual(rcalc.get_rate_BA(), 1955395816., delta=1e4)
        
        rcalc = RateCalculation(self.db.transition_states(), [m1,m3], [m2, m4], T=0.592)
        rcalc.compute_rates()
        self.assertAlmostEqual(rcalc.get_rate_AB(), 8638736600., delta=1e4)
        self.assertAlmostEqual(rcalc.get_rate_BA(), 3499625167., delta=1e4)
        
        rla = RatesLinalg(self.db.transition_states(), [m1], [m2], T=0.592)
        rAB = rla.compute_rates()
        self.assertAlmostEqual(rAB, 7106337458., delta=1e4)
        
        rla = RatesLinalg(self.db.transition_states(), [m1, m3], [m2, m4], T=0.592)
        rAB = rla.compute_rates()
        self.assertAlmostEqual(rAB, 8638736600., delta=1e4)
        


if __name__ == "__main__":
    unittest.main()

