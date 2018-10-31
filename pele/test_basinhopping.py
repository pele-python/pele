import unittest
from numpy import abs

from pele.basinhopping import BasinHopping
from pele.systems import LJCluster

class TestBasinhopping(unittest.TestCase):
    def setUp(self):
        natoms = 6
        self.system = LJCluster(natoms)
        self.e1 = -12.7120622568
        self.e2 = -12.302927529580728
    
    def assertEnergy(self, e):
        self.assertTrue( abs(e-self.e1) < 1e-4 or abs(e-self.e2) < 1e-4)
    
    def test_create_basinhopping(self):
        bh = self.system.get_basinhopping(outstream=None)
        bh.run(3)
        self.assertEnergy(bh.result.energy)
        self.assertEnergy(bh.markovE)

    def test_takestep(self):
        takestep = self.system.get_takestep(stepsize=1.)
        bh = self.system.get_basinhopping(outstream=None, takestep=takestep)
        bh.run(3)
        self.assertEnergy(bh.result.energy)
        self.assertEnergy(bh.markovE)
        
        
if __name__ == "__main__":
    unittest.main()
