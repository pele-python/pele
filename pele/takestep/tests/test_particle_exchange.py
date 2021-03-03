import unittest
import numpy as np

from pele.takestep import ParticleExchange


class TestParticleExchange(unittest.TestCase):
    def setUp(self):
        Alist = [1, 3, 5]
        Blist = [0, 2, 4]
        x = np.zeros([6, 3], int)
        for i in Alist:
            x[i, :] = [1, 1, 1]

        self.step = ParticleExchange(Alist, Blist)
        self.x = x.flatten()

    def test(self):
        x = self.x
        self.step(x)
        x = x.reshape(-1, 3)
        sA = x[self.step.Alist, :].sum()
        self.assertEqual(sA, 6)
        sB = x[self.step.Blist, :].sum()
        self.assertEqual(sB, 3)


if __name__ == "__main__":
    unittest.main()
        
