import numpy as np
import unittest

class TestMinDist(unittest.TestCase):
    """
    a base class for mindist unit tests
    """
    def runtest(self, X1, X2, mindist):
        X1i = np.copy(X1)
        X2i = np.copy(X2)
        
        (distreturned, X1, X2) = mindist(X1, X2)

        distinit = np.linalg.norm(self.X1 - X2i)
        distfinal = np.linalg.norm(X1 - X2)
        self.assertAlmostEqual( distfinal, distreturned, 5, "returned distance is wrong: %g != %g, difference %g" % (distfinal, distreturned, distfinal-distreturned) )
        self.assertLessEqual(distfinal, distinit)
        
        #test if the energies have changed
        Ei = self.pot.getEnergy(X1i)        
        Ef = self.pot.getEnergy(X1)
        self.assertAlmostEqual( Ei, Ef, 10, "Energy of X1 changed: %g - %g = %g" % (Ei, Ef, Ei - Ef) )
        Ei = self.pot.getEnergy(X2i)        
        Ef = self.pot.getEnergy(X2)
        self.assertAlmostEqual( Ei, Ef, 10, "Energy of X2 changed: %g - %g = %g" % (Ei, Ef, Ei - Ef) )

        return distreturned, X1, X2

    
    

