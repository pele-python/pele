from __future__ import print_function
from __future__ import absolute_import
import numpy as np

from pele.potentials import BasePotential
from .fortran import lj as ljf


__all__ = ["LJ"]


class LJ(BasePotential):
    """ simple lennard jones potential"""

    def __init__(self, eps=1.0, sig=1.0, boxl=None):
        self.sig = sig
        self.eps = eps
        self.boxl = boxl
        if self.boxl is None:
            self.periodic = False
            self.boxl = 10000.
        else:
            self.periodic = True

    def getEnergy(self, coords):
        E = ljf.ljenergy(
            coords, self.eps, self.sig, self.periodic, self.boxl)
        return E

    def getEnergyGradient(self, coords):
        E, grad = ljf.ljenergy_gradient(
            coords, self.eps, self.sig, self.periodic, self.boxl)
        return E, grad

    def getEnergyList(self, coords, ilist):
        E = ljf.energy_ilist(
            coords, self.eps, self.sig, ilist.reshape(-1), self.periodic,
            self.boxl)
        return E

    def getEnergyGradientList(self, coords, ilist):
        E, grad = ljf.energy_gradient_ilist(
            coords, self.eps, self.sig, ilist.reshape(-1), self.periodic,
            self.boxl)
        return E, grad

    def getEnergyGradientHessian(self, coords):
        if self.periodic: raise Exception("Hessian not implemented for periodic boundaries")
        from .fortran.lj_hess import ljdiff

        g, energy, hess = ljdiff(coords, True, True)
        return energy, g, hess


#
# only testing below here
#

def main():  # pragma: no cover
    # test class
    natoms = 12
    coords = np.random.uniform(-1, 1, natoms * 3) * 2

    lj = LJ()
    E = lj.getEnergy(coords)
    print("E", E)
    E, V = lj.getEnergyGradient(coords)
    print("E", E)
    print("V")
    print(V)

    print("try a quench")
    from pele.optimize import mylbfgs as quench

    quench(coords, lj, iprint=1)
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )


if __name__ == "__main__":
    main()

