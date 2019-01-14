from __future__ import print_function
from __future__ import absolute_import
import unittest

import numpy as np

from pele.landscape import DoubleEndedConnect
from .test_graph import create_random_database
from pele.systems import LJCluster


class TestDistanceGraph(unittest.TestCase):
    def setUp(self):
        nmin = 10
        natoms=13
        
        sys = LJCluster(natoms)
        
        pot = sys.get_potential()
        mindist = sys.get_mindist()
        
        db = create_random_database(nmin=nmin, natoms=natoms, nts=nmin//2)
        min1, min2 = list(db.minima())[:2]
        
        
        self.db = db
        self.natoms = natoms
        self.connect = DoubleEndedConnect(min1, min2, pot, mindist, db,
                                          merge_minima=True, max_dist_merge=1e100)
#        for ts in db.transition_states():
#            self.add_ts(ts.energy, ts.coords, 
#                        ts.minimum1.energy, ts.minimum1.coords,
#                        ts.minimum2.energy, ts.minimum2.coords)
        for m in db.minima():
            self.connect.dist_graph.addMinimum(m)
    
    def make_result(self, coords, energy):
        from pele.optimize import Result
        res = Result()
        res.coords = coords
        res.energy = energy
        res.eigenval = 1.
        res.eigenvec = coords.copy()
        return res

    def add_ts(self, tse, tsx, m1e, m1x, m2e, m2x):
        tsres = self.make_result(tsx, tse)
        min_ret1 = self.make_result(m1x, m1e)
        min_ret2 = self.make_result(m2x, m2e)
        return self.connect._addTransitionState(tsres, min_ret1, min_ret2)

    
    def test_merge_minima(self):
        # merge two minima and make sure the distance graph is still ok
        min3, min4 = list(self.db.minima())[2:4]
        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "the distance graph is broken at the start")

        print(min3.id(), min4.id(), "are connected", self.connect.graph.areConnected(min3, min4))
        print(min3.id(), "number of edges", self.connect.graph.graph.degree(min3))
        print(min4.id(), "number of edges", self.connect.graph.graph.degree(min4))
        self.connect.mergeMinima(min3, min4)
        
        self.assertNotIn(min4, self.connect.graph.graph)
        self.assertNotIn(min4, self.connect.dist_graph.Gdist)
        self.assertNotIn(min4, self.db.minima())
        
        allok = self.connect.dist_graph.checkGraph()
        
        
        self.assertTrue(allok, "merging broke the distance graph")
        
    def test_add_TS_existing_minima(self):
        min3, min4 = list(self.db.minima())[4:6]
        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "the distance graph is broken at the start")

        print(min3.id(), min4.id(), "are connected", self.connect.graph.areConnected(min3, min4))
        print(min3.id(), "number of edges", self.connect.graph.graph.degree(min3))
        print(min4.id(), "number of edges", self.connect.graph.graph.degree(min4))
        
        coords = np.random.uniform(-1,1,self.natoms*3)
        E = float(min3.energy + min4.energy)
        tsres = self.make_result(coords, E)
        min_ret1 = self.make_result(min3.coords, min3.energy)
        min_ret2 = self.make_result(min4.coords, min4.energy)
        
       
        self.connect._addTransitionState(tsres, min_ret1, min_ret2)

        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "adding a transition state broke the distance graph")
        

    def test_add_TS_new_minima(self):
        min3 = list(self.db.minima())[6]
        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "the distance graph is broken at the start")

#        print min3.id(), min4.id(), "are connected", self.connect.graph.areConnected(min3, min4)
        print(min3.id(), "number of edges", self.connect.graph.graph.degree(min3))

        #create new minimum from thin air
        coords = np.random.uniform(-1,1,self.natoms*3)
        E = np.random.rand()*10.
        min_ret1 = self.make_result(min3.coords, min3.energy)
        min_ret2 = self.make_result(coords, E)

#        min_ret2 = [coords, E]
#        min_ret1 = [min3.coords, min3.energy]


        #create new TS from thin air        
        coords = np.random.uniform(-1,1,self.natoms*3)
        E = float(min3.energy + min_ret2.energy)
        tsres = self.make_result(coords, E)
        
        self.connect._addTransitionState(tsres, min_ret1, min_ret2)

        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "adding a transition state broke the distance graph")

    def run_add_TS(self, min3, min4, nocheck=False):
        if not nocheck:
            allok = self.connect.dist_graph.checkGraph()
            self.assertTrue(allok, "the distance graph is broken at the start")

        print(min3.id(), min4.id(), "are connected", self.connect.graph.areConnected(min3, min4))
        print(min3.id(), "number of edges", self.connect.graph.graph.degree(min3))
        print(min4.id(), "number of edges", self.connect.graph.graph.degree(min4))
        
        coords = np.random.uniform(-1,1,self.natoms*3)
        E = float(min3.energy + min4.energy)
        tsres = self.make_result(coords, E)
        min_ret1 = self.make_result(min3.coords, min3.energy)
        min_ret2 = self.make_result(min4.coords, min4.energy)

#        min_ret1 = [min3.coords, min3.energy]
#        min_ret2 = [min4.coords, min4.energy]
        
        self.connect._addTransitionState(tsres, min_ret1, min_ret2)

        if not nocheck:
            allok = self.connect.dist_graph.checkGraph()
            self.assertTrue(allok, "adding a transition state broke the distance graph")


    def test_add_TS_existing_not_connected(self):
        minima = list(self.db.minima())
        min3 = minima[2]
        for min4 in minima[3:]:
            if not self.connect.graph.areConnected(min3, min4):
                break
        self.run_add_TS(min3, min4)
        
    def test_add_TS_existing_already_connected(self):
        minima = list(self.db.minima())
        min3 = minima[2]
        for min4 in minima[3:]:
            if self.connect.graph.areConnected(min3, min4):
                break
        self.run_add_TS(min3, min4)
    
    def test_add_TS_multiple(self):
        minima = list(self.db.minima())
        min3 = minima[2]
        nnew = 4
        for min4 in minima[3:3+nnew]:
            self.run_add_TS(min3, min4, nocheck=True)

        allok = self.connect.dist_graph.checkGraph()
        self.assertTrue(allok, "adding multiple transition states broke the distance graph")

if __name__ == "__main__":
    unittest.main()

