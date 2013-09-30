import unittest

from pele.systems import LJCluster
from pele.transition_states import GeneralizedDimer
from pele.utils.xyz import read_xyz

class TestGeneralizedDimer(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(18)
        self.pot = self.system.get_potential()
    
    def make_dimer(self, x):
        return GeneralizedDimer(x, self.pot,
                                leig_kwargs=dict(orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors()
                                                 ) ,
                                )
    
    def test1(self):
        x = self.system.get_random_configuration()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)

    def test2(self):
        xyz = read_xyz(open("lj18_ts.xyz", "r"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)
#        print res

if __name__ == "__main__":
    unittest.main()
