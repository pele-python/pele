from __future__ import print_function
from __future__ import absolute_import
import numpy as np

from pele.potentials import BasePotential
from pele.systems import BaseSystem


class Booth(BasePotential):
    target_E = 0.
    target_coords = np.array([1., 3.])
    xmin = np.array([-10., -10.])
    # xmin = np.array([0., 0.])
    xmax = np.array([10., 10.])

    def getEnergy(self, coords):
        x, y = coords
        return (x + 2.*y - 7.)**2 + (2.*x + y - 5.)**2


class BoothSystem(BaseSystem):
    def get_potential(self):
        return Booth()

    def get_random_configuration(self, eps=1e-3):
        pot = self.get_potential()
        xmin, xmax = pot.xmin, pot.xmax
        x = np.random.uniform(xmin[0] + eps, xmax[0] - eps)
        y = np.random.uniform(xmin[1] + eps, xmax[1] - eps)
        return np.array([x, y])


def test1():
    from ._beale import makeplot2d

    s = BoothSystem()
    f = s.get_potential()
    f.test_potential(f.target_coords)
    print("")
    f.test_potential(s.get_random_configuration())
    f.test_potential(np.array([1., 1.]))  # , print_grads=True)

    # from base_function import makeplot2d
    makeplot2d(f, nx=60, zlim=[0, 100])


if __name__ == "__main__":
    test1()

