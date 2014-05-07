import unittest
import os
import nose

import numpy as np
from pele.mindist import PointGroupOrderCluster, ExactMatchAtomicCluster

class TestPgorderLj75(unittest.TestCase):
    """as of Mar 5 2014 this test fails.  It needs to be fixed"""
    def test1(self):
        d = os.path.dirname(__file__)
        fname = os.path.join(d, "coords.lj75.gmin")
        coords = np.genfromtxt(fname)
        print fname
        self.assertEqual(coords.size, 75*3)
        
        permlist = [range(75)]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(coords)
        print pgorder
        
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
        print m.energy
        permlist = [range(6)]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(m.coords)
        print pgorder
        
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
        permlist = [range(natoms)]
        match = ExactMatchAtomicCluster(permlist=permlist, can_invert=True)
        calculator = PointGroupOrderCluster(match)
        pgorder = calculator(m.coords)
        print pgorder
        
        self.assertEqual(pgorder, 120)

        
        
    
if __name__ == "__main__":
    unittest.main() 
