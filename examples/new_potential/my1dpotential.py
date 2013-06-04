"""
a 1d example potential
"""
import numpy as np

from pele.potentials import BasePotential

class My1DPot(BasePotential):
    """1d potential"""
    def getEnergy(self, x):
        return np.cos(14.5 * x[0] - 0.3) + (x[0] + 0.2) * x[0]

    def getEnergyGradient(self, x):
        E = self.getEnergy(x)
        grad = np.array([-14.5 * np.sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2])
        return E, grad

#
#system class
#

from pele.systems import BaseSystem
class My1DSystem(BaseSystem):
    def get_potential(self):
        return My1DPot()


def run_basinhopping():
    sys = My1DSystem()
    database = sys.create_database()
    x0 = np.array([1.])
    print x0[:]
    bh = sys.get_basinhopping(database=database, coords=x0)
    bh.run(100)
    print "found", len(database.minima()), "minima"
    min0 = database.minima()[0]
    print "lowest minimum found at", min0.coords, "with energy", min0.energy


if __name__ == "__main__":
    import numpy as np
    pot = My1DPot()
    coords = np.array([1.])
    e = pot.getEnergy(coords)
    print e
    e, g = pot.getEnergyGradient(coords)
    print e
    
    gnum = pot.NumericalDerivative(coords, eps=1e-6)
    print np.max(np.abs(gnum-g)), np.max(np.abs(gnum))
    print np.max(np.abs(gnum-g)) / np.max(np.abs(gnum))

    run_basinhopping()