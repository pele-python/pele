from __future__ import absolute_import
import unittest
import numpy as np
import os

from pele.potentials._lj_cpp import BLJCut
from pele.utils.xyz import read_xyz
from pele.utils.rotations import vec_random_ndim
from pele.potentials.ljpshiftfast import LJpshift
from . import _base_test

_x0 = np.array([ 2.77463041, -2.08121793, -0.04984229,  0.95964575,  2.6318472 ,
       -2.74566884, -3.85139038, -3.43982285, -0.10923911,  0.85063569,
        0.5508115 , -1.46110073,  3.90892924,  0.63796175, -0.95887062,
        0.40758575,  1.96267545,  1.35386315, -1.88064354, -3.46932132,
       -1.03932642,  1.03774006, -2.31860792,  2.02204443, -3.46770815,
       -1.91747921,  2.43803651, -2.45252574,  1.11568705,  0.19736247,
        3.39846376, -1.89362584, -3.47231127,  1.88052771,  2.17742424,
        3.26252682,  3.45577655, -3.88838742, -2.12510331,  0.93422686,
        3.59213057,  3.60140895,  0.45322551,  3.3248508 ,  1.13252967,
       -0.87993829, -0.11207466,  0.83448386,  0.39638337,  3.40945141,
        3.34986749, -0.8409951 ,  3.70610023, -2.60835467, -2.98936384,
       -2.91936674,  0.04529733, -3.82780156,  3.58376169,  2.61692377,
       -3.87984815, -2.59042996, -1.34349141, -2.95202524,  2.47592554,
       -1.24210678,  3.52085986,  0.65611344,  3.03065588,  2.75787556,
        3.24313855, -0.32095787,  0.37077453,  2.38882873, -1.71424919,
       -0.07797182,  0.79288246, -3.8757338 ,  0.74785127, -0.53058921,
        2.45888423, -1.47804158,  3.14310967,  0.62285772, -2.52791839,
        2.30343387,  0.89624942, -3.56872582, -0.63845056,  1.43255069,
        3.34881422, -3.9967838 ,  3.81407319, -0.98735748,  3.79026831,
        0.83772881,  2.63076646,  0.59769204,  1.02460959, -1.71538975,
        0.69466673,  2.00017411,  2.86651069,  2.04065751,  1.58445799,
        2.91583544, -1.41855203,  1.36631033, -0.39300851, -0.94317798,
       -0.7135092 , -0.78816333, -1.46092843,  0.97535494, -0.55802183,
        3.79041662,  1.42240713, -2.41144089, -0.58639193, -1.25323008,
        2.38111043,  3.03998631,  3.23073565,  1.3017585 , -1.8383339 ,
       -1.98106639,  2.83918354,  0.22171717,  2.41728867,  0.57990814,
        1.8651402 ,  0.15209302,  2.16707128,  0.55086393, -0.27432097,
       -1.25848874, -3.45432521, -0.97660657, -3.36299138,  3.86253691,
       -2.54709719,  2.49486958,  2.99969316,  1.50730602,  0.5559553 ,
       -2.71222851, -0.26495982, -1.23862359, -2.19968034,  0.74009495])

class TestBLJ_CPP_simple(_base_test._TestConfiguration):
    def setUp(self):
        natoms = 50
        rcut = 1.6
        ntypeA = int(natoms * 0.8)
        self.x0 = _x0.copy()
        self.pot = BLJCut(natoms, ntypeA, rcut=rcut)

        self.e0 = 1412.0144910476681
        
        self.pot_fortran = LJpshift(natoms, ntypeA, rcut=rcut)
    
    def test_against_fortran(self):
        efort = self.pot_fortran.getEnergy(self.x0)
        self.assertAlmostEqual(self.e0, efort)

class TestBLJ_CPP_Bulk(_base_test._TestConfiguration):
    def setUp(self):
        natoms = 50
        rcut = 1.6
        ntypeA = int(natoms * 0.8)
        self.x0 = _x0.copy()
        boxl = 7.
        boxvec = np.array([boxl]*3)
        self.pot = BLJCut(natoms, ntypeA, rcut=rcut, boxvec=boxvec)
        self.pot_fortran = LJpshift(natoms, ntypeA, boxl=boxl, rcut=rcut)
        self.e0 = 6971.685336750815
    
    def test_against_fortran(self):
        efort = self.pot_fortran.getEnergy(self.x0)
        self.assertAlmostEqual(self.e0, efort)

class TestBLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        np.random.seed(1)
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_blj13_min.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        ntypeA, self.Emin, rcut, epsAA, sigAA, epsBB, sigBB, epsAB, sigAB = list(map(float, xyz.title.split()[1::2]))
        ntypeA = int(ntypeA)
        self.rcut = rcut

        natoms = self.xmin.size // 3

        self.pot = BLJCut(natoms, ntypeA, rcut=rcut, sigAA=sigAA, epsAA=epsAA,
                          epsBB=epsBB, sigBB=sigBB, epsAB=epsAB, sigAB=sigAB)
        self.xrandom = np.random.uniform(-1, 1, self.xmin.size) * 5.


class TestBLJ_CPP_Cut(unittest.TestCase):
    def setUp(self):
        np.random.seed(1)
        self.rcut = 2.5
        self.atom1 = np.random.uniform(-1, 1, [3])
        v = vec_random_ndim(3)
        v /= np.linalg.norm(v)
        self.v = v

    def test_rcut(self):
        pot = BLJCut(2, 1, rcut=self.rcut)
        atom2 = self.atom1 + self.rcut * 1.001 * self.v
        x = np.array(list(self.atom1) + list(atom2))
        e = pot.getEnergy(x)
        self.assertLess(e, 1e-20)

    def test_rcut2(self):
        pot = BLJCut(2, 1, rcut=self.rcut)
        atom2 = self.atom1 + self.rcut * 0.7 * self.v
        # print np.linalg.norm(atom2 - self.atom1)
        x = np.array(list(self.atom1) + list(atom2))
        e = pot.getEnergy(x)
        self.assertGreater(np.abs(e), 1e-5)


def makeplot():  # pragma: no cover
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
    rlist = np.arange(.85, 3, .01)
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
    # makeplot()
    unittest.main()

