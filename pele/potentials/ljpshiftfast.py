from __future__ import absolute_import
from .fortran import ljpshiftfort as ljpshiftfort

from pele.potentials import BasePotential

__all__ = ["LJpshift"]


class BLJ_interaction_type:
    """
    holds the parameters for a given interaction type: AA, AB, BB
    """

    def __init__(self, eps, sig, rcut):
        self.eps = eps
        self.sig = sig
        self.rcut = rcut * self.sig

        self.ircut2 = 1.0 / self.rcut ** 2
        self.sig6 = self.sig ** 6
        self.sig12 = self.sig6 ** 2

        sigrc6 = self.sig6 / self.rcut ** 6
        sigrc12 = sigrc6 ** 2

        self.const = 4.0 * sigrc6 - 7.0 * sigrc12
        self.rconst = (6.0 * sigrc12 - 3.0 * sigrc6) / self.rcut ** 2


class LJpshift(BasePotential):
    """binary lennard jones potential with smooth cutoff"""

    def __init__(self, natoms, ntypeA, boxl=None, rcut=2.5, epsBB=0.5, sigBB=0.88, epsAB=1.5, sigAB=0.8):
        self.boxl = boxl
        self.natoms = natoms
        self.ntypeA = ntypeA

        sigAA = 1.
        epsAA = 1.
        self.AA = BLJ_interaction_type(epsAA, sigAA, rcut)
        self.BB = BLJ_interaction_type(epsBB, sigBB, rcut)
        self.AB = BLJ_interaction_type(epsAB, sigAB, rcut)

        if self.boxl is None:
            self.periodic = False
            self.boxl = 10000.
        else:
            self.periodic = True

    def getEnergy(self, coords):
        V, E = ljpshiftfort.ljpshift(coords, False, False,
                                     self.boxl, self.boxl, self.boxl,
                                     self.AA.rcut, self.periodic, self.ntypeA,
                                     self.AB.eps, self.BB.eps, self.AB.sig, self.BB.sig)
        return E

    def getEnergyGradient(self, coords):
        V, E = ljpshiftfort.ljpshift(coords, True, False,
                                     self.boxl, self.boxl, self.boxl,
                                     self.AA.rcut, self.periodic, self.ntypeA,
                                     self.AB.eps, self.BB.eps, self.AB.sig, self.BB.sig)
        return E, V

