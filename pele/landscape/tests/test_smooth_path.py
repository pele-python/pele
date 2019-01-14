import unittest


import numpy as np

from pele.landscape import smooth_path

def crap_mindist(x1, x2):
    dist = np.linalg.norm(x1 - x2)
    return dist, x1, x2

class TestSmoothPath(unittest.TestCase):
    def test1(self):
        natoms = 13
        m1 = np.random.uniform(-1,1,natoms*3)
        ts = np.random.uniform(-1,1,natoms*3)
        m2 = np.random.uniform(-1,1,natoms*3)
        path = [m1, ts, m2]
        
        density = 5.
        spath = smooth_path(path, crap_mindist, density=density)
        
        d1 = crap_mindist(m1,ts)[0]
        d2 = crap_mindist(ts,m2)[0]
        dmax = max(d1, d2)
        for x1, x2 in zip(spath[:-1], spath[1:]):
            d = crap_mindist(x1, x2)[0]
            self.assertGreater(d, 0)
            self.assertLess(d, dmax)
        
        self.assertGreaterEqual(len(spath), 2*density)
        
        
         

if __name__ == "__main__":
    unittest.main()

