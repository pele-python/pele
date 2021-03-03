import unittest

from pele.takestep import AdaptiveStepsize, AdaptiveStepsizeTemperature, RandomDisplacement
from pele.systems import LJCluster


class TestTakestepAdaptiveStepsize(unittest.TestCase):
    def test1(self):
        system = LJCluster(6)
        ss0 = 1.
        displace = RandomDisplacement(stepsize=ss0)
        ts = AdaptiveStepsize(displace, interval=10)
        bh = system.get_basinhopping(takestep=ts)
        bh.run(9)
        self.assertAlmostEqual(displace.stepsize, ss0, 10)
        bh.run(2)
        self.assertNotAlmostEqual(displace.stepsize, ss0, 1)


class TestTakestepAdaptiveStepTemperature(unittest.TestCase):
    def test1(self):
        system = LJCluster(6)
        ss0 = 1.
        displace = RandomDisplacement(stepsize=ss0)
        ts = AdaptiveStepsizeTemperature(displace, interval=10)
        bh = system.get_basinhopping(takestep=ts)
        bh.run(9)
        self.assertAlmostEqual(displace.stepsize, ss0, 10)
        bh.run(2)
        self.assertNotAlmostEqual(displace.stepsize, ss0, 1)


if __name__ == "__main__":
    unittest.main()

