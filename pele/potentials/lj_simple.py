from __future__ import print_function
import numpy as np

from pele.potentials import BasePotential


class LJ(BasePotential):
    """
    A stripped down, simple version of LJ for example purposes
    """

    def __init__(self, eps=1.0, sig=1.0):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps

    def vij(self, r):
        return 4. * self.eps * ( (self.sig / r) ** 12 - (self.sig / r) ** 6 )

    def dvij(self, r):
        return 4. * self.eps * ( -12. / self.sig * (self.sig / r) ** 13 + 6. / self.sig * (self.sig / r) ** 7 )

    def getEnergy(self, coords):
        natoms = coords.size / 3
        coords = np.reshape(coords, [natoms, 3])
        energy = 0.
        for i in range(natoms):
            for j in range(i + 1, natoms):
                dr = coords[j, :] - coords[i, :]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
        return energy

    def getEnergyGradient(self, coords):
        natoms = coords.size / 3
        coords = np.reshape(coords, [natoms, 3])
        energy = 0.
        V = np.zeros([natoms, 3])
        for i in range(natoms):
            for j in range(i + 1, natoms):
                dr = coords[j, :] - coords[i, :]
                r = np.linalg.norm(dr)
                energy += self.vij(r)
                g = self.dvij(r)
                V[i, :] += -g * dr / r
                V[j, :] += g * dr / r
        V = V.reshape([natoms * 3])
        return energy, V


#
# testing only below here
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

    ret = quench(coords, lj, iprint=-1)
    #quench( coords, lj.getEnergyGradientNumerical, iprint=1 )
    print("energy ", ret.energy)
    print("rms gradient", ret.rms)
    print("number of function calls", ret.nfev)


if __name__ == "__main__":
    main()

