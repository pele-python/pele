from __future__ import print_function
import unittest
import os
import sys
import nose

import numpy as np
from pele.mindist import PointGroupOrderCluster, ExactMatchAtomicCluster
from pele.utils.xyz import read_xyz


class TestPgorderLj75(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        d = os.path.dirname(__file__)
        fname = os.path.join(d, "coords.lj75.gmin.xyz")
        xyz = read_xyz(open(fname, "r"))
        coords = xyz.coords.reshape(-1)
        print(fname)
        self.assertEqual(coords.size, 75*3)
        
        permlist = [list(range(75))]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(coords)
#        print pgorder
        
        self.assertEqual(pgorder, 20)

class TestPgorderLj6(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        from pele.systems import LJCluster
        system = LJCluster(6)
        db = system.create_database()
        bh = system.get_basinhopping(db)
        bh.setPrinting(ostream=None)
        while db.minima()[0].energy > -12.7:
            bh.run(10)

        m = db.minima()[0]
        print(m.energy)
        permlist = [list(range(6))]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(m.coords)
#        print pgorder
        
        self.assertEqual(pgorder, 48)

class TestPgorderLj13(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        from pele.systems import LJCluster
        natoms = 13
        system = LJCluster(natoms)
        db = system.create_database()
        bh = system.get_basinhopping(db)
        bh.setPrinting(ostream=None)
        while db.minima()[0].energy > -44.3:
            bh.run(10)

        m = db.minima()[0]
        permlist = [list(range(natoms))]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(m.coords)
#        print pgorder
        
        self.assertEqual(pgorder, 120)

class TestPgorderLj13Database(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        d = os.path.dirname(__file__)
        dbfname = os.path.join(d, "lj13_small_pathsample.{}.sqlite".format(sys.version_info.major))

        from pele.systems import LJCluster
        natoms = 13
        system = LJCluster(natoms)
        db = system.create_database(dbfname, createdb=False)

        permlist = [list(range(natoms))]
        
        ts_min = list(db.minima()) + list(db.transition_states())
        
        for m in ts_min:
            match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
            calculator = PointGroupOrderCluster(match)
            pgorder = calculator(m.coords)
            self.assertEqual(pgorder, m.pgorder)
#            print pgorder
        

class TestPgorderLj75Database(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        d = os.path.dirname(__file__)
        dbfname = os.path.join(d, "lj75_very_small_pathsample.{}.sqlite".format(sys.version_info.major))

        from pele.systems import LJCluster
        natoms = 75
        system = LJCluster(natoms)
        db = system.create_database(dbfname, createdb=False)

        permlist = [list(range(natoms))]
        
        ts_min = list(db.minima()) + list(db.transition_states())
        
        for m in ts_min:
            match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
            calculator = PointGroupOrderCluster(match)
            pgorder = calculator(m.coords)
            self.assertEqual(pgorder, m.pgorder)
#            print pgorder
        

        
    
if __name__ == "__main__":
    unittest.main() 

