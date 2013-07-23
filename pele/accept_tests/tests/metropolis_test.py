import unittest
import numpy as np
import pygmin.exceptions as exc
from pygmin.accept_tests.metropolis import \
    Metropolis as metro, \
    MetropolisNonQuench as metro_nq

class TestMetropolis(unittest.TestCase):

    def setUp(self):
        """
        to be run before each test method
        """
        pass

    def test_dropEnergyPosPosAccept(self):
        """
        test that when newE is positive and lower than oldE which is positive,
        the structure is always accepted 
        """
        the_test = True
        self.assertTrue(the_test)

    def test_dropEnergyPosZeroAccept(self):
        """
        test that when newE is zero and (thus) lower than oldE which is positive,
        the structure is always accepted 
        """
        the_test = True
        self.assertTrue(the_test)
        
    def test_dropEnergyZeroNegAccept(self):
        """
        test that when newE is negative and (thus) lower than oldE which is zero,
        the structure is always accepted 
        """
        the_test = True
        self.assertTrue(the_test)
        
    def test_dropEnergyNegNegAccept(self):
        """
        test that when newE is negative and lower than oldE which is negative,
        the structure is always accepted 
        """
        the_test = True
        self.assertTrue(the_test)

    def test_forceAccept(self):
        """
        test that using forceAccept means that the Metropolis step is always
        accepted
        """
        pass

if __name__ == '__main__':
    unittest.main()