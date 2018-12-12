from __future__ import print_function
import unittest
from copy import deepcopy

import numpy as np

from pele.angleaxis.molecules import create_water
from pele.angleaxis import RBTopology, TransformAngleAxisCluster, MeasureAngleAxisCluster
import pele.angleaxis.aamindist as am
from pele.utils.rotations import random_aa, aa2mx
from pele.utils import rotations
from pele.potentials.tests._base_test import assert_arrays_almost_equal
from pele.angleaxis.aatopology import AASiteType
from pele.angleaxis.rigidbody import RigidFragment

class TestAAMindist(unittest.TestCase):
    def setUp(self):
        #GMIN.initialize()
#        self.pot = GMINPotential(GMIN)
        self.nrigid = 10
        self.water = create_water()
        self.topology = RBTopology()
        self.topology.add_sites([deepcopy(self.water) for _ in range(self.nrigid)])
        self.topology.finalize_setup()

#    def test_zeroev(self):
#        x = self.pot.getCoords()
#        zev = self.topology.zeroEV(x)
#        
#        eps = 1e-5
#        for dx in zev:    
#            print "ev test", (self.pot.getEnergy(x) - self.pot.getEnergy(x + eps*dx))/eps  
#        
#        dx = np.random.random(x.shape)
#        dx/=np.linalg.norm(dx)
#        print "ev test", (self.pot.getEnergy(x) - self.pot.getEnergy(x + eps*dx))/eps
        
    def test_distance(self):  
        for i in range(100):
            coords1 = np.random.random(6*self.nrigid)*4
            coords2 = np.random.random(6*self.nrigid)*4
            
            coords1[:3*self.nrigid]=0
            coords2[:3*self.nrigid]=0
            
            measure1 = am.MeasureAngleAxisCluster(self.topology)
            measure2 = am.MeasureRigidBodyCluster(self.topology)
            
            #print 
            self.assertAlmostEqual(measure1.get_dist(coords1, coords2),
                                   measure2.get_dist(coords1, coords2))
    

class TestAATransform(unittest.TestCase):
    def setUp(self):
        self.nrigid = 10
        self.topology = RBTopology()
        self.topology.add_sites([create_water() for _ in range(self.nrigid)])
        self.topology.finalize_setup()
        self.transform = TransformAngleAxisCluster(self.topology)
    
    def test_cpp_rotate(self):
        x0 = np.array([random_aa() for _ in range(2*self.nrigid)]).ravel()
        aa = random_aa()
        mx = aa2mx(aa)
        
        x0_rotated = x0.copy()
        self.transform.rotate(x0_rotated, mx)
        
        x0_rotated_python = x0.copy()
        self.transform._rotate_python(x0_rotated_python, mx)
        
        assert_arrays_almost_equal(self, x0_rotated, x0_rotated_python)
        
        mx_reverse = aa2mx(-aa)
        self.transform.rotate(x0_rotated, mx_reverse)
        assert_arrays_almost_equal(self, x0, x0_rotated)

def create_tetrahedron():
    f = RigidFragment()
    f.add_atom("h", [ 1., 0., -1./np.sqrt(2)])
    f.add_atom("h", [-1., 0., -1./np.sqrt(2)])
    f.add_atom("h", [0.,  1., 1./np.sqrt(2)])
    f.add_atom("h", [0., -1., 1./np.sqrt(2)])
    f.finalize_setup()
    return f

class TestAAMeasure(unittest.TestCase):
    def setUp(self):
        self.nrigid = 10
        self.topology = RBTopology()
        self.topology.add_sites([create_tetrahedron() for _ in range(self.nrigid)])
        self.topology.finalize_setup()
        self.measure = MeasureAngleAxisCluster(self.topology)
    
    def test_cpp_align(self):
        x1 = np.array([random_aa() for _ in range(2*self.nrigid)]).ravel()
        x2 = np.array([random_aa() for _ in range(2*self.nrigid)]).ravel()
        
        x1_cpp = x1.copy()
        x2_cpp = x2.copy()
        x1_python = x1.copy()
        x2_python = x2.copy()
        
        ret = self.measure._align_pythonic(x1_python, x2_python)
        self.assertIsNone(ret)
        ret = self.measure.align(x1_cpp, x2_cpp)
        self.assertIsNone(ret)
        
        assert_arrays_almost_equal(self, x1_python, x1)
        assert_arrays_almost_equal(self, x1_cpp, x1)
        
        # assert that the center of mass coordinates didn't change
        assert_arrays_almost_equal(self, x2_python[:3*self.nrigid], x2[:3*self.nrigid])
        assert_arrays_almost_equal(self, x2_cpp[:3*self.nrigid], x2[:3*self.nrigid])
        
        # assert that x2 actually changed
        max_change = np.max(np.abs(x2_cpp - x2))
        self.assertGreater(max_change, 1e-3)

        assert_arrays_almost_equal(self, x2_python, x2_cpp)

    def test_align_bad_input(self):
        x1 = np.array([random_aa() for _ in range(2*self.nrigid)]).ravel()
        x2 = list(x1)
        
        with self.assertRaises(TypeError):
            self.measure.align(x1, x2)

    def test_symmetries(self):
        tet = create_tetrahedron()
        print(tet.symmetries)
        self.assertEqual(len(tet.symmetries), 12)


    def test_align_exact(self):
        x1 = np.array([random_aa() for _ in range(2*self.nrigid)]).ravel()
        x2 = x1.copy()
        tet = create_tetrahedron()
        mx = tet.symmetries[2].copy() 
        p = x2[-3:].copy()
        x2[-3:] = rotations.rotate_aa(rotations.mx2aa(mx), p)
        
        self.measure._align_pythonic(x1, x2)
        assert_arrays_almost_equal(self, x1, x2)
        




    
if __name__ == '__main__':
    unittest.main()

