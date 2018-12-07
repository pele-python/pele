from __future__ import print_function
from __future__ import absolute_import
import numpy as np

from pele.potentials import LJ
from pele.potentials import BasePotential
from .fortran import AT as ATfort

__all__ = ["ATLJ"]


class ATLJ(BasePotential):
    """
    Lennard Jones + three body Axilrod-Teller term
    
    V = sum_ij VLJ_ij   +  sum_ijk  Z * (1 + 3*cos(t1)*cos(t2)*cos(t3)) / (rij * rjk * rik)**3 )
    
    where t1, t2, t3 are the internal angles of the triangle ijk
    
    Z > 0 stabilizes linear vs. triangular geometries 
    """

    def __init__(self, eps=1.0, sig=1.0, Z=1.):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps
        self.Z = Z
        self.lj = LJ(self.sig, self.eps)

    def getEnergySlow(self, coords):
        Elj = self.lj.getEnergy(coords)

        natoms = coords.size // 3
        X = np.reshape(coords, [natoms, 3])
        Z = self.Z
        energy = 0.
        for i in range(natoms):
            for j in range(i):
                for k in range(j):
                    drij = X[i, :] - X[j, :]
                    drik = X[i, :] - X[k, :]
                    drjk = X[j, :] - X[k, :]
                    rij = np.linalg.norm(drij)
                    rik = np.linalg.norm(drik)
                    rjk = np.linalg.norm(drjk)
                    energy += (Z * (1. + 3. *
                                    np.dot(drij, -drjk) *
                                    np.dot(-drij, -drik) *
                                    np.dot(drik, drjk) / (rij * rik * rjk) ** 2)
                               / (rij * rik * rjk) ** 3 )
        energy += Elj
        return energy

    def getEnergyFortran(self, coords):
        garbage, e = ATfort.axt(coords, False, self.Z)
        Elj = self.lj.getEnergy(coords)
        return e + Elj

    def getEnergyGradientFortran(self, coords):
        grad, e = ATfort.axt(coords, True, self.Z)

        elj, gradlj = self.lj.getEnergyGradient(coords)
        return e + elj, grad + gradlj

    def getEnergy(self, coords):
        return self.getEnergyFortran(coords)

    def getEnergyGradient(self, coords):
        return self.getEnergyGradientFortran(coords)


#
# testing only below here
#

def testing():  # pragma: no cover
    # test class
    natoms = 3
    coords = np.random.uniform(-1, 1, natoms * 3) * 2

    lj = ATLJ(Z=1.)

    E = lj.getEnergy(coords)
    print("E", E)
    E, V = lj.getEnergyGradient(coords)
    print("E", E)
    print("V")
    print(V)

    print("try a quench")
    from pele.optimize import mylbfgs as quench

    ret = quench(coords, lj, iprint=-1)
    print("energy ", ret.energy)
    print("rms gradient", ret.rms)
    print("number of function calls", ret.nfev)

    from pele.utils.xyz import write_xyz

    printlist = []
    for i in range(100):
        coords = np.random.uniform(-1, 1, natoms * 3) * 2
        ret = quench(coords, lj.getEnergyGradient, iprint=-1)
        coords = ret.coords
        X = np.reshape(coords, [natoms, 3])
        com = X.sum(0) / natoms
        X[:, :] -= com[np.newaxis, :]
        printlist.append(np.reshape(X, natoms * 3))

    with open("out.xyz", "w") as fout:
        for coords in printlist:
            write_xyz(fout, coords)


if __name__ == "__main__":
    testing()

