import unittest
import os

from pele.systems import LJCluster
from pele.transition_states._generalized_dimer import GeneralizedDimer
from pele.utils.xyz import read_xyz


class TestGeneralizedDimer(unittest.TestCase):
    def setUp(self):
        self.system = LJCluster(18)
        self.pot = self.system.get_potential()

    def make_dimer(self, x):
        return GeneralizedDimer(
            x,
            self.pot,
            leig_kwargs=dict(
                orthogZeroEigs=self.system.get_orthogonalize_to_zero_eigenvectors()
            ),
        )

    # TODO: this fails for ~0.15% of random starts (3 of 2000; np.random.seed 10089, 10735
    # and 11584 before get_random_configuration reproduce it). In each, an atom breaks off
    # the cluster and the dimer converges to a minimum plus a free atom (curvature ~ +1e-8),
    # since the dimer has no guard against dissociation. Also, res.success only checks that
    # the lowest curvature found is negative: 56 of the same 2000 runs "succeed" at points
    # that the exact Hessian shows are not first-order saddles (mostly two negative modes).
    def test1(self):
        x = self.system.get_random_configuration()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)

    def test2(self):
        # get the path of the file directory
        path = os.path.dirname(os.path.abspath(__file__))
        xyz = read_xyz(open(path + "/lj18_ts.xyz"))
        x = xyz.coords.flatten()
        dimer = self.make_dimer(x)
        res = dimer.run()
        self.assertTrue(res.success)
        self.assertLess(res.nsteps, 5)


#        print res

if __name__ == "__main__":
    unittest.main()
