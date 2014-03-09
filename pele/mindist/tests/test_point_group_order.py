import unittest
import os
import nose

import numpy as np
from pele.mindist import PointGroupOrderCluster, ExactMatchAtomicCluster

@nose.tools.nottest # this test is know to fail.  remove this line when it passes
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

        
        
    
if __name__ == "__main__":
    unittest.main() 
