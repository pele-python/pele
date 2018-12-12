import unittest


import numpy as np

from pele.mindist.permutational_alignment import  find_permutations_munkres, \
    find_permutations_OPTIM, find_best_permutation, find_permutations_hungarian


class PermutationTest(unittest.TestCase):
    #    def __init__(self, *args):
    coords = np.array([[0., 0., 1.], [0., 0., 2.], [0., 0, 3.]])
    
    def test_identity(self):
        self.check_perm([0, 1, 2])
        
    def test_simple(self):
        self.check_perm([0, 2, 1])
    
    def test_non_invertable(self):
        self.check_perm([1, 2, 0])
    
    def check_perm(self, perm):
        coords1 = self.coords.copy()
        coords2 = self.coords.copy()[perm]
        dist, perm_calc = find_best_permutation(coords1, 
                              coords2,
                              reshape=False,
                              user_algorithm = find_permutations_OPTIM)
        try:
            self.assertItemsEqual(perm, perm_calc)
        except AttributeError:
            self.assertCountEqual(perm, perm_calc)
        self.assertAlmostEqual(np.linalg.norm(coords1 - self.coords), 0.)
        self.assertAlmostEqual(np.linalg.norm(coords2[perm_calc] - self.coords), 0.)
    
    def test_munkres(self):
        perm = [1, 0, 2]
        coords1 = self.coords.copy()
        coords2 = self.coords.copy()[perm]
        dist1, perm_calc1 = find_best_permutation(coords1.copy(), 
                              coords2,
                              reshape=False,
                              user_algorithm=find_permutations_hungarian)
        dist2, perm_calc2 = find_best_permutation(coords1.copy(), 
                              coords2,
                              reshape=False,
                              user_algorithm=find_permutations_munkres)
        self.assertEqual(perm_calc1, perm_calc2)
        
class PermutationTestBinary(unittest.TestCase):
    #    def __init__(self, *args):
    coords = np.array([[0., 0., 1.], [0., 0., 2.], [0., 0, 3.], [0., 0, 4.]])
    permlist = [[0, 2], [1, 3]]
       
    def test_identity(self):
        self.check_perm([0, 1, 2, 3])
    def test_perm(self):
        self.check_perm([2, 1, 0, 3])
        self.check_perm([2, 3, 0, 1])
        self.check_perm([0, 3, 2, 1])

    def test_nonpermutable(self):
        coords1 = self.coords.copy()
        coords2 = self.coords.copy()[[0, 2, 1, 3]]
        dist, perm_calc = find_best_permutation(coords1, 
                              coords2,
                              reshape=False,
                              permlist=self.permlist,
                              user_algorithm = find_permutations_OPTIM)
        try:
            self.assertItemsEqual(perm_calc, [0, 2, 1, 3])
        except AttributeError:
            self.assertCountEqual(perm_calc, [0, 2, 1, 3])
                
    def check_perm(self, perm):
        coords1 = self.coords.copy()
        coords2 = self.coords.copy()[perm]
        dist, perm_calc = find_best_permutation(coords1, 
                              coords2,
                              reshape=False,
                              permlist=self.permlist,
                              user_algorithm = find_permutations_OPTIM)
        try:
            self.assertItemsEqual(perm, perm_calc)
        except AttributeError:
            self.assertCountEqual(perm, perm_calc)
        self.assertAlmostEqual(np.linalg.norm(coords1 - self.coords), 0.)
        self.assertAlmostEqual(np.linalg.norm(coords2[perm] - self.coords), 0.)

        
if __name__ == "__main__":
    unittest.main()

