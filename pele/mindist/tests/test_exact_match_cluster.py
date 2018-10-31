import unittest
import numpy as np
import itertools

from pele.mindist import ExactMatchAtomicCluster
from pele.systems import LJCluster
import pele.utils.rotations as rotations

class TestExactMatchAtomicCluster(unittest.TestCase):
    def setUp(self):
        self.natoms = 11
        self.system = LJCluster(self.natoms)
        self.match = ExactMatchAtomicCluster(permlist=self.system.get_permlist(), can_invert=True)
    
    def rc(self):
        return self.system.get_random_configuration()
    
    def rtrans(self, x):
        vec3 = np.random.uniform(-1,1,3)
        self.match.transform.translate(x, vec3)
    
    def rrot(self, x):
        aa = rotations.random_aa()
        mx = rotations.aa2mx(aa)
        self.match.transform.rotate(x, mx)
    
    def invert(self, x):
        self.match.transform.invert(x)
        
    
    def rperm(self, x):
        perm = list(range(self.natoms))
        np.random.shuffle(perm)
        self.match.transform.permute(x, perm)
        
    
    def are_same_check(self, tform):
        x1 = self.rc()
        x2 = x1.copy()
        
        tform(x1)
        tform(x2)
        
        x1_bkup = x1.copy()
        x2_bkup = x2.copy()
        
        self.assertTrue(self.match(x1, x2))
        self.assertTrue((x1==x1_bkup).all())
        self.assertTrue((x2==x2_bkup).all())
        
#        transform = self.match.find_transformation(x1, x2)
#        self.assertIsNotNone(transform, "the structures should compare exact")
#        self.assertTrue((x1==x1_bkup).all())
#        self.assertTrue((x2==x2_bkup).all())
#        x2_tformed = x2.copy()
#        self.match.apply_transformation(x2_tformed, transform)
#        self.assertAlmostEqual(np.linalg.norm(x1_bkup - x2_tformed), 0., 3)
        
    def find_transform_check(self, tform):
        x1 = self.rc()
        x2 = x1.copy()
        
        tform(x1)
        tform(x2)
        
        x1_bkup = x1.copy()
        x2_bkup = x2.copy()
        
        transform = self.match.find_transformation(x1, x2)
#        print "invert", transform.invert
#        print "translate", transform.translation
#        print "permute", transform.permutation
#        print "rotate", transform.rotation
        self.assertIsNotNone(transform, "the structures should compare exact")
        self.assertTrue((x1==x1_bkup).all(), "the passed structures should not be modified")
        self.assertTrue((x2==x2_bkup).all(), "the passed structures should not be modified")
        # check that transform is the transformation that turns x2 into x1
        x2_tformed = x2.copy()
        self.match.apply_transformation(x2_tformed, transform)
        self.assertAlmostEqual(np.linalg.norm(x1_bkup - x2_tformed), 0., 3)

    def apply_checks(self, tform):
        self.are_same_check(tform)
        self.find_transform_check(tform)

    def test_translate(self):
        self.apply_checks(self.rtrans)

    def test_rotate(self):
        self.apply_checks(self.rrot)
        
    def test_permute(self):
        self.apply_checks(self.rperm)
        
    def test_invert(self):
        self.apply_checks(self.rperm)
    
    def fchain(self, flist):
        """return a function which applies all of the functions in flist to the input"""
        def function_chain(x):
            for f in reversed(flist):
                f(x)
        return function_chain
    
    def test_2(self):
        # run apply_checks() on all combinations of length 2 of the transformations
        flist_all = [self.invert, self.rrot, self.rtrans, self.rperm]
        for flist in itertools.product(flist_all, repeat=2):
            tform = self.fchain(flist)
            self.apply_checks(tform)
        
    def test_3(self):
        # run apply_checks() on all combinations of lenth 3 of the transformations
        flist_all = [self.invert, self.rrot, self.rtrans, self.rperm]
        for flist in itertools.product(flist_all, repeat=3):
            tform = self.fchain(flist)
            self.apply_checks(tform)
        
    def test_4(self):
        # run apply_checks() on all combinations of lenth 4 of the transformations
        flist_all = [self.invert, self.rrot, self.rtrans, self.rperm]
        for flist in itertools.product(flist_all, repeat=4):
            tform = self.fchain(flist)
            self.apply_checks(tform)
        
        

if __name__ == "__main__":
    unittest.main()
        

