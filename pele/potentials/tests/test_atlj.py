from past.utils import old_div
import unittest
import numpy as np
from pele.potentials import LJ, ATLJ


class TestATLJ(unittest.TestCase):
    def testenergy(self):
        natoms = 10
        coords = np.random.uniform(-1, 1, natoms * 3) * 2

        from pele.optimize import mylbfgs as quench

        lj = LJ()
        ret = quench(coords, lj)
        coords = ret.coords

        atlj = ATLJ(Z=3.0)
        e2 = atlj.getEnergySlow(coords)
        # e1 = atlj.getEnergyWeave(coords)
        # #print "%g - %g = %g" % (e1, e2, e1-e2)
        # self.assertTrue( abs(e1 - e2) < 1e-12, "ATLJ: two energy methods give different results: %g - %g = %g" % (e1, e2, e1-e2) )

        e1 = atlj.getEnergyFortran(coords)
        # print "%g - %g = %g" % (e1, e2, e1-e2)
        # print e1/e2
        self.assertTrue(
            abs(e1 - e2) < 1e-12,
            "ATLJ: fortran energy gives different results: %g - %g = %g"
            % (e1, e2, e1 - e2),
        )

    def testGradient(self):
        natoms = 10
        coords = np.random.uniform(-1, 1, natoms * 3) * 2

        atlj = ATLJ(Z=3.0)

        e, Gf = atlj.getEnergyGradientFortran(coords)
        Gn = atlj.NumericalDerivative(coords)
        print(Gf)
        print(Gn)
        maxdiff = np.max(np.abs(Gf - Gn))
        maxnorm = old_div(np.max(np.abs(Gf + Gn)), 2)
        maxrel = np.max(np.abs(old_div((Gf - Gn), (Gf + Gn) * 2.0)))
        print(
            "maximum relative difference in gradients",
            maxdiff,
            old_div(maxdiff, maxnorm),
        )
        self.assertTrue(
            old_div(maxdiff, maxnorm) < 1e-4,
            "ATLJ: gradient differs from numerical gradient by %g" % maxdiff,
        )


if __name__ == "__main__":
    unittest.main()
