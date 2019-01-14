from __future__ import print_function
import numpy as np

from pele.potentials import BasePotential
from pele.potentials.fortran import scdiff_periodic as fortran_sc


class SuttonChen(BasePotential):
    """The sutton chen potential
    
    for modelling the surfaces of metal crystals
    
    :
        First calculate the potential energy. Choosing SIG=sqrt(2) makes the
        unit of length the bulk nearest-neighbour distance. Choosing SIG to 
        be the bulk lattice constant (=sqrt(2)x nn-distance) gives the original
        SC form. Can also choose sigma=eps=1 and then scale the energies by eps
        and the distances by sig.
        metal  NN  MM   SCEPS/eV    SCC   SIG
         Cu     9   6   0.012382  39.432
         Ni     9   6   0.015707  39.432
         Ag    12   6   0.002542  144.41
         Au    10   8   0.012793  34.408
         Pt    10   8   0.019833  34.408
         Pd    12   7   0.004179  108.27
    

    """

    def __init__(self, eps=1., sig=1., c=144.41, boxvec="default", rcut=2.5,
                 n=10, m=8):
        if boxvec == "default":
            boxvec = [10., 10., 10.]
        if boxvec is None:
            raise ValueError("non-periodic sutton chen potential not implemented")
        self.eps = eps
        self.sig = sig
        self.c = c
        self.n = n
        self.m = m
        self.boxvec = np.array(boxvec, dtype=float)
        self.rcut = rcut

    def getEnergy(self, coords):
        e, g = self.getEnergyGradient(coords)
        return e

    def getEnergyGradient(self, coords):
        g, e = fortran_sc.scdiff_periodic(coords, self.eps, self.c, self.sig,
                                          self.n, self.m, self.boxvec, self.rcut)
        return e, g


#
# testing only below here
#

def test_silver():  # pragma: no cover
    rcut = 2.5
    potnocut = SuttonChen(rcut=30000.2, boxvec=[10000] * 3, c=144.41, n=12, m=6)
    pot = SuttonChen(rcut=rcut, boxvec=[10000] * 3, c=144.41, n=12, m=6)
    natoms = 2
    x = np.random.uniform(-2, 2, [3 * natoms])
    pot.getEnergy(x)

    pot.test_potential(x)

    import matplotlib.pyplot as plt

    x0 = .5
    x = np.array([0.] * 5 + [x0])
    pot.test_potential(x)
    energies = []
    energies_nocut = []
    dxlist = np.arange(0, 5, .0001)
    for dx in dxlist:
        x1 = x.copy()
        x1[-1] += dx
        energies.append(pot.getEnergy(x1))
        energies_nocut.append(potnocut.getEnergy(x1))
    plt.plot(dxlist + x0, energies, label="rcut=%f" % rcut)
    plt.plot(dxlist + x0, energies_nocut, label="no cutoff")
    plt.legend()
    plt.show()


def test_platinum():  # pragma: no cover
    rcut = 2.5
    pot = SuttonChen(rcut=rcut, boxvec=[10000] * 3, c=34.408, n=10, m=8)
    natoms = 20
    x = np.random.uniform(-2, 2, [3 * natoms])
    pot.getEnergy(x)
    pot.test_potential(x)


def test_fcc():  # pragma: no cover
    from pele.utils.xyz import read_xyz

    fname = "fcc100-coords.8x8x8.xyz"
    # pot = SuttonChen(rcut=3.2, boxvec=[8.] * 3, c=144.41, n=12, m=6)
    pot = SuttonChen(rcut=3.2, boxvec=[8.] * 3, c=34.408, n=10, m=8)
    xyz = read_xyz(open(fname, "r"))
    x = xyz.coords.flatten()

    e, g = pot.getEnergyGradient(x)
    print("energy", e)
    print("norm grad", np.linalg.norm(g))
    print("rms grad", np.linalg.norm(g) / np.sqrt(g.size))


if __name__ == "__main__":
    test_fcc()
# test_silver()
# test_platinum()

