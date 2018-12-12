from __future__ import division
from __future__ import absolute_import
import unittest
import os
import logging

import numpy as np

from pele.potentials import _hs_wca_cpp
from pele.optimize._quench import lbfgs_cpp
from . import _base_test


def read_xyzdr(fname, bdim=3):
    coords = []
    radii = []
    rattlers = []
    f = open(fname, "r")
    while True:
        xyzdr = f.readline()
        if not xyzdr: break
        x, y, z, d, r = xyzdr.split()
        coords.extend([float(x), float(y), float(z)])
        radii.extend([float(d) / 2])
        for _ in range(bdim):
            rattlers.extend([float(r)])
    return np.array(coords), np.array(radii), np.array(rattlers)


def minimize(coords, pot):
    result = lbfgs_cpp(coords, pot)
    # result = modifiedfire_cpp(coords, pot)
    return result.coords, result.energy, result.grad, result.rms


class TestHS_WCA_CPP(_base_test._BaseTest):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        xyz, hs_radii, rattlers = read_xyzdr(current_dir + "/_hswca20_min2.xyzdr")
        sca = 0.205071132088
        boxv = [6.26533756282, 6.26533756282, 6.26533756282]
        self.pot = _hs_wca_cpp.HS_WCA(eps=1, sca=sca, radii=hs_radii, boxvec=boxv, use_periodic=True)
        self.natoms = 20
        result = minimize(xyz, self.pot)
        self.xmin = result[0]  # xyz
        self.Emin = result[1]  # self.pot.getEnergy(self.xmin)
        self.xrandom = np.random.uniform(-1, 1, len(xyz)) * 1e-2

# class TestHS_WCA_CPP_NeighborList(_base_test._BaseTest):
# def setUp(self):
#        current_dir = os.path.dirname(__file__)
#        xyz, hs_radii, rattlers = read_xyzdr(current_dir + "/_hswca20_min.xyzdr")
#        sca = 0.259921049895
#        boxv= [6.85206773233, 6.85206773233, 6.85206773233]
#        self.Emin = 182.943079825
#        np.random.seed(0)
#        self.natoms = 20
#        nlist = [[i,j] for i in xrange(self.natoms) for j in xrange(i+1,self.natoms)]
#        nlist = np.array(nlist, dtype=np.int64).reshape(-1)
#        self.pot = _hs_wca_cpp.HS_WCANeighborList(nlist, eps=1, sca=sca, radii=hs_radii)
#        self.xrandom = np.random.uniform(-1,1,len(xyz)) *1e-2
#        self.xmin = xyz


if __name__ == "__main__":
    logging.basicConfig(filename='hs_wca_cpp.log', level=logging.DEBUG)
    unittest.main()

