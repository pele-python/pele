from __future__ import print_function
from __future__ import absolute_import
import unittest
import numpy as np
import os
import logging

from pele.potentials import _lj_cpp
from pele.utils.xyz import read_xyz
from . import _base_test


class TestLJ_CPP(_base_test._BaseTest):
    def setUp(self):
        self.pot = _lj_cpp.LJ()
        self.natoms = 13
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)


class TestErrorPotential(unittest.TestCase):
    def setUp(self):
        self.pot = _lj_cpp._ErrorPotential()
        self.x = np.random.uniform(-1, 1, [9])

    def test(self):
        with self.assertRaises(RuntimeError):
            self.pot.getEnergy(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.getEnergyGradient(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalDerivative(self.x)
        with self.assertRaises(RuntimeError):
            self.pot.NumericalHessian(self.x)


# with self.assertRaises(NotImplementedError):
# pot.getEnergyGradient(_xrand)

class TestLJ_CPP_NeighborList(_base_test._BaseTest):
    def setUp(self):
        np.random.seed(0)
        self.natoms = 13
        nlist = [[i, j] for i in range(self.natoms) for j in range(i + 1, self.natoms)]
        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
        self.pot = _lj_cpp.LJNeighborList(nlist)
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_lj13_gmin.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin = float(xyz.title)


class TestLJCutCellLists(_base_test._TestConfiguration):
    def setUp(self):
        boxv = np.array([10.] * 3)
        natoms = 30
        rcut = 2.
        self.x0 = np.array([ 1.39757936,  1.83281398,  2.32259516,  1.60560537,  3.11214842,
                            2.61582972,  4.76275273,  3.37651599,  4.66292907,  6.12541448,
                            3.75706966,  0.54669172,  4.55566734,  7.41983351,  3.0333311 ,
                            3.5405557 ,  0.71801575,  1.71847253,  2.17861941,  5.70588003,
                            6.74642534,  3.43588101,  2.0687562 ,  5.55499349,  3.27452166,
                            3.96526201,  6.61476763,  0.64229015,  1.27522777,  1.21875002,
                            0.99191894,  4.41664435,  1.97658992,  6.41180268,  2.15398194,
                            4.2636531 ,  6.95210635,  2.75332174,  7.29174054,  6.53110874,
                            4.64505199,  6.90914585,  2.9151178 ,  6.28565612,  7.95207857,
                            6.06415512,  4.92514773,  2.53285987,  3.72928997,  0.53255714,
                            6.51884144,  1.2042502 ,  7.34367396,  4.25806453,  2.06642627,
                            4.83650925,  1.29926411,  3.12454566,  2.71078146,  2.99731328])
        
        # use a random configuration that is quenched slightly so the energy is not crazy
        x = np.random.uniform(0, boxv[0], natoms*3)
        self.pot_true = _lj_cpp.LJCut(boxvec=boxv, rcut=rcut)
        self.pot = _lj_cpp.LJCutCellLists(boxvec=boxv, rcut=rcut, ncellx_scale=1.)
        from pele.optimize import lbfgs_cpp
        ret = lbfgs_cpp(x, self.pot_true, tol=10.)
        self.x0 = ret.coords
        
#         print repr(self.x0)
        self.e0 = self.pot_true.getEnergy(self.x0)
        print("true energy", self.e0)

        self.ae_kwargs = dict(places=6)
        
        self.rcut = rcut
        self.boxvec = boxv
    
    def check_cell_density(self, ncellx_scale=1.):
        # note: 11/2014: this becomes extremely slow for ncellx_scale > 2.
        self.pot = _lj_cpp.LJCutCellLists(boxvec=self.boxvec, rcut=self.rcut, ncellx_scale=ncellx_scale)
        
        e = self.pot.getEnergy(self.x0)
        
        self.assertAlmostEqual(e, self.e0, **self.ae_kwargs)
        
        
    def test_cell_density(self):
        for f in [0.1, 0.5, 1., 2.]:
            self.check_cell_density(ncellx_scale=f)


if __name__ == "__main__":
    logging.basicConfig(filename='lj_cpp.log', level=logging.DEBUG)
    unittest.main()

