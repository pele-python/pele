import unittest

from pele.optimize import _quench
from pele.systems import LJCluster


class _PotWrapper(object):
    def __init__(self, pot):
        self.nfev = 0
        self.pot = pot

    def getEnergy(self, coords):
        self.nfev += 1
        return self.pot.getEnergy(coords)

    def getGradient(self, coords):
        self.nfev += 1
        return self.pot.getGradient(coords)

    def getEnergyGradient(self, coords):
        self.nfev += 1
        return self.pot.getEnergyGradient(coords)


class TestMinimizers_NFEV(unittest.TestCase):
    def setUp(self):
        natoms = 31
        self.system = LJCluster(natoms)
        self.pot = _PotWrapper(self.system.get_potential())

        # get a partially minimized structure
        self.x0 = self.system.get_random_configuration()

        self.minimizers = [_quench.lbfgs_py,
                           _quench.mylbfgs,
                           _quench.lbfgs_scipy,
        ]


    def do_check(self, minimizer):
        self.pot.nfev = 0
        ret = minimizer(self.x0, self.pot)
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertGreater(ret.nfev, 0)


    def test1(self):
        self.do_check(_quench.lbfgs_py)

    def test2(self):
        self.do_check(_quench.mylbfgs)

    def test3(self):
        self.do_check(_quench.fire)

    def test4(self):
        self.do_check(_quench.lbfgs_scipy)

    def test5(self):
        self.do_check(_quench.bfgs_scipy)

    def test6(self):
        self.do_check(_quench.cg)


if __name__ == "__main__":
    unittest.main()
        