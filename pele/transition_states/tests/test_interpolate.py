import unittest

import numpy as np

from pele.mindist.periodic_exact_match import MeasurePeriodic
from pele.transition_states._interpolate import InterpolateLinearMeasure


class TestInterpolatePeriodic(unittest.TestCase):
    def test1(self):
        natoms = 10
        L = 4.
        boxvec = np.ones(3) * L
        xi = np.random.uniform(0,L,3*natoms)
        xf = np.random.uniform(0,L,3*natoms)
        
        measure = MeasurePeriodic(boxvec)
        
        dist = measure.get_dist(xi, xf)
        
        int = InterpolateLinearMeasure(measure)
        
        for f in [0, .1, .2, .3, .8, .9, 1]:
            x3 = int(xi, xf, f)
            d2 = measure.get_dist(xi, x3)
            self.assertAlmostEqual(f*dist, d2)

if __name__ == "__main__":
    unittest.main()
        
