import unittest
import numpy as np
import os

from pele.potentials._lj_cpp import BLJCut
from pele.utils.xyz import read_xyz
from pele.utils.rotations import vec_random_ndim
import _base_test

class TestBLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_blj13_min.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        ntypeA, self.Emin, rcut, epsAA, sigAA, epsBB, sigBB, epsAB, sigAB = map(float, xyz.title.split()[1::2])
        ntypeA = int(ntypeA)
        self.rcut = rcut
        
        natoms = self.xmin.size / 3
        
        self.pot = BLJCut(natoms, ntypeA, rcut=rcut, sigAA=sigAA, epsAA=epsAA, 
                          epsBB=epsBB, sigBB=sigBB, epsAB=epsAB, sigAB=sigAB) 
        self.xrandom = np.random.uniform(-1,1,self.xmin.size) *5.

class TestBLJ_CPP_Cut(unittest.TestCase):
    def setUp(self):
        self.rcut = 2.5
        self.atom1 = np.random.uniform(-1,1,[3])
        v = vec_random_ndim(3)
        v /= np.linalg.norm(v)
        self.v = v 

    def test_rcut(self):
        pot = BLJCut(2, 1, rcut=self.rcut)
        atom2 = self.atom1 + self.rcut * 1.001  * self.v
        x = np.array(list(self.atom1) + list(atom2))
        e = pot.getEnergy(x)
        self.assertLess(e, 1e-20)

    def test_rcut2(self):
        pot = BLJCut(2, 1, rcut=self.rcut)
        atom2 = self.atom1 + self.rcut * 0.7 * self.v
#        print np.linalg.norm(atom2 - self.atom1)
        x = np.array(list(self.atom1) + list(atom2))
        e = pot.getEnergy(x)
        self.assertGreater(np.abs(e), 1e-5)
        
        
def makeplot():
    atom1 = np.zeros(3)
    rcut = 2.5
    potAA = BLJCut(2, 0, rcut=rcut)
    potAB = BLJCut(2, 1, rcut=rcut)
    potBB = BLJCut(2, 2, rcut=rcut)
    eAA = []
    eAB = []
    eBB = []
    v = vec_random_ndim(3)
    v /= np.linalg.norm(v)
    rlist = np.arange(.85,3,.01)
    for r in rlist:
        atom2 = atom1 + r * v
        x = np.array(list(atom1) + list(atom2))
        eAA.append(potAA.getEnergy(x))
        eAB.append(potAB.getEnergy(x))
        eBB.append(potBB.getEnergy(x))

    import matplotlib.pyplot as plt
    plt.plot(rlist, eAA)
    plt.plot(rlist, eAB)
    plt.plot(rlist, eBB)
    plt.show()
    

if __name__ == "__main__":
#    makeplot()
    unittest.main()
