"""
Created on 6 Apr 2012

@author: ruehle
"""
from __future__ import print_function

import numpy as np
import copy

__all__ = ["QuenchBenchmark"]


class PotentialWrapper(object):
    def __init__(self, potential):
        self.potential = potential
        self.reset()

    def reset(self):
        self.energies = []

    def getEnergy(self, coords):
        E = self.potential.getEnergy(coords)
        self.energies.append(E)
        return E

    def getEnergyGradient(self, coords):
        E, g = self.potential.getEnergyGradient(coords)
        self.energies.append(E)
        return E, g

    def getGradient(self, coords):
        E, g = self.potential.getEnergyGradient(coords)
        self.energies.append(E)
        return g


class QuenchBenchmark(object):
    """
    classdocs
    """


    def __init__(self, potential):
        """
        Constructor
        """
        self.potential = PotentialWrapper(potential)
        self.minimizer = []

    def addMinimizer(self, label, minimizer):
        self.minimizer.append([label, minimizer, 0.0, None])

    def run(self, Emin, coords):
        for minimizer in self.minimizer:
            self.potential.reset()
            print("Testing Minimizer " + minimizer[0])

            E, grad = self.potential.getEnergyGradient(coords)
            res = minimizer[1](coords, self.potential)
            minimizer[2] = res.energy
            minimizer[3] = np.array(self.potential.energies).copy() - Emin
            print("Minimizer " + minimizer[0] + ": " + str(res.energy))

    def plot(self):
        import pylab as pl

        for m in self.minimizer:
            pl.loglog(np.array(m[3]), label=m[0])

        pl.legend()
        pl.xlabel("energy evaluations")
        pl.ylabel("energy")
        pl.show()


if __name__ == "__main__":
    import pele.potentials.lj as lj
    import scipy.optimize
    from pele.optimize import _quench as quench

    print("Running benchmark with lennard jones potential")
    pot = lj.LJ()

    natoms = 36

    coords = np.random.random(3 * natoms) * 10.0
    res = quench.lbfgs_scipy(coords, pot, tol=1e-3)
    coords = res.coords
    coords += np.random.random(coords.shape) * 0.1
    res = quench.lbfgs_scipy(coords, pot, tol=1e-3)
    Emin = res.energy

    bench = QuenchBenchmark(pot)
    bench.addMinimizer("lbfgs", quench.lbfgs_scipy)
    bench.addMinimizer("mylbfgs", quench.mylbfgs)
    bench.addMinimizer("lbfgs_py", quench.lbfgs_py)
    bench.addMinimizer("cg", quench.cg)
    bench.addMinimizer("fire", quench.fire)
    bench.addMinimizer("bfgs_scipy", quench.bfgs_scipy)

    print("The reference energy is " + str(Emin))
    bench.run(Emin, coords)
    bench.plot()

