import unittest

import numpy as np

from pele.optimize import _quench
from pele.systems import LJCluster

_x0 = np.array([-0.31257026,  0.93852964, -0.75764933,  1.10735569, -0.28162525,
               -0.93779627, -1.13033058,  0.71269579,  1.46762412,  2.04395519,
                1.33642043,  0.92190506, -0.64272146, -1.88548667, -0.43575506,
                1.53658785,  0.50666837, -0.78503166,  1.12270223, -0.75936947,
                0.63422829, -1.46169159, -0.90870614,  2.15510392, -1.83216771,
               -1.97126243, -2.07119939,  0.36629137,  2.26810183, -0.35254211,
               -2.18638657,  0.19291511,  0.22061239,  0.75544507,  1.42411619,
               -0.34456975,  1.86810091, -1.46365119,  0.55227977,  1.75721169,
               -0.93030996, -0.97132648,  1.24388712,  0.04651615,  0.10942041,
                1.8919052 ,  0.60505881,  0.40692323,  0.64136677, -0.7701783 ,
               -0.3544892 ,  2.20478373, -0.24784546,  0.18697057, -1.5483192 ,
                0.04081523, -1.10145187, -0.6199117 ,  0.37408135, -1.43531892,
               -1.8423539 ,  1.81237162,  1.10705189,  0.95055242,  1.37036181,
               -2.18625576,  2.12650874,  1.77119623, -0.45147687, -0.26442857,
               -0.30932327,  0.0631677 , -1.91862341,  1.46364406, -0.19226111,
                1.23098819, -1.58045323, -1.36607109, -0.49081348,  1.82354102,
               -1.57057896, -1.47315425, -1.01904413, -0.86166884, -2.15676312,
                2.18839054, -1.59410816, -1.82811369,  2.18126752, -0.78893266,
               -1.45326636, -1.13734021, -1.94326349])

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
        self.x0 = _x0.copy()

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
        # note: on rare occasions this raises the error:
        # glibc detected *** /usr/bin/python: double free or corruption
        self.do_check(_quench.bfgs_scipy)

    def test6(self):
        self.do_check(_quench.cg)


if __name__ == "__main__":
    unittest.main()
        
