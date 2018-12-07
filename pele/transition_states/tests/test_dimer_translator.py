import unittest
import numpy as np

from pele.systems import LJCluster
from pele.transition_states._dimer_translator import _DimerTranslator
from pele.utils.rotations import vec_random_ndim




class TestTransverseWalker_NFEV(unittest.TestCase):
    def setUp(self):
        from pele.optimize.tests.test_nfev import _PotWrapper
        self.system = LJCluster(18)
        self.pot = _PotWrapper(self.system.get_potential())
    
    def test(self):
        x = self.system.get_random_configuration()
        evec = vec_random_ndim(x.size)
        evec /= np.linalg.norm(evec) 
        opt = _DimerTranslator(x, self.pot, evec)
        ret = opt.run(10)
        self.assertEqual(ret.nfev, self.pot.nfev)
        self.assertGreater(ret.nfev, 0)
 


if __name__ == "__main__":
    unittest.main()

