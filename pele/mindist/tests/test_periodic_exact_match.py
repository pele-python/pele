import unittest
import copy
import numpy as np

from pele.mindist.periodic_exact_match import MeasurePeriodic, ExactMatchPeriodic

class TestExactMatchPeriodicLJ(unittest.TestCase):
    def setUp(self):
        self.natoms = 100
        rho = .5
        boxl = (float(self.natoms) / rho)**(1./3)
        boxlengths = np.ones(3) * boxl + np.random.rand(3)*.1
        
        self.permlist = self.get_permlist()        
        self.measure = MeasurePeriodic(boxlengths, self.permlist)
        
        self.x1 = self.get_random_configuration()
        self.x2diff = self.get_random_configuration()
        self.x2same = self.randomly_permute(self.x1.copy())
        
        self.exact_match = ExactMatchPeriodic(self.measure, accuracy=1e-5)
    
    def get_permlist(self):
        return [list(range(self.natoms))]
    
    def randomly_permute(self, x):
        import random
        x = x.reshape(-1,3)
        xnew = x.copy()
        for atomlist in self.permlist:
            permutation = copy.copy(atomlist)
            random.shuffle(permutation)
            xnew[atomlist,:] = x[permutation,:]
        return xnew.flatten()

    
    def get_random_configuration(self):
        return np.random.uniform(-1,1,self.natoms*3) * (float(self.natoms))**(1./3) * 1.5
    
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.measure.boxlengths  
        self.assertTrue(self.exact_match(self.x1, self.x2same))
    
    def test_exact_match_translate(self):
        dx = np.random.uniform(-1,1,3)
        self.x2same = (self.x1.reshape(-1,3) + dx[np.newaxis,:]).ravel()
        self.assertTrue(self.exact_match(self.x1, self.x2same))

    def test_exact_match_translate_permute(self):
        dx = np.random.uniform(-1,1,3)
        self.x2same = (self.x1.reshape(-1,3) + dx[np.newaxis,:]).ravel()
        self.x2same = self.randomly_permute(self.x2same)
        self.assertTrue(self.exact_match(self.x1, self.x2same))

    def test_exact_match_permute_translate(self):
        dx = np.random.uniform(-1,1,3)
        self.x2same = self.randomly_permute(self.x1.copy())
        self.x2same = (self.x2same.reshape(-1,3) + dx[np.newaxis,:]).ravel()
        self.assertTrue(self.exact_match(self.x1, self.x2same))


class TestExactMatchPeriodicBLJ(TestExactMatchPeriodicLJ):        
    def get_permlist(self):
        self.ntypeA = int(self.natoms *.8)
        return [list(range(self.ntypeA)), list(range(self.ntypeA, self.natoms))]
    
if __name__ == '__main__':
    unittest.main()

