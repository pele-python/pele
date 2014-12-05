import unittest
from itertools import izip

import numpy as np
from numpy import sqrt, cos, sin, pi

from pele.angleaxis import RBTopology, RigidFragment, RBPotentialWrapper
from pele.potentials import LJ
from pele.angleaxis._otp_bulk import OTPBulk
from pele.thermodynamics import get_thermodynamic_information
from pele.utils import rotations
from pele.angleaxis._aa_utils import _rot_mat_derivative, _sitedist_grad, _sitedist
from pele.angleaxis.aamindist import MeasureRigidBodyCluster

def put_in_box(x, boxvec):
    x = x.reshape(-1, boxvec.size)
    x -= boxvec * np.round(x / boxvec)

_x1 = np.array([13.248024173529, -2.518155646687, -0.690721764444, 
                -5.416744538918, 6.710532681251, -4.497894715783, 
                9.711353031813, -7.301769746043, -20.474698913868, 
                -6.295179676919, 1.928520894483, 14.934303628184, 
                -15.548248597546, 1.213746414430, 4.229236560578, 
                0.783370927992, -2.404788961421, -0.890713258205, 
                3.419810220016, -0.071815622421, -0.760530821079, 
                -1.143135910803, 1.554601397812, 0.519659901046, 
                2.070362891228, -2.206974019735, -1.059841980313, 
                0.366427346535, 2.352197191462, 3.355031578983])

_x2 = np.array([-2.006975474436, -9.197065954327, -6.590822695704, 
                1.800970988133, 6.368229090117, -1.404656537018, 
                1.702333632220, -8.977566443712, 2.737913663382, 
                -3.213668111055, 10.316330897574, -0.704781226877, 
                -2.583456515026, 1.522946980579, -0.536961031572, 
                -1.303060588871, 1.809537546480, -1.108327987283, 
                -0.369995110332, 0.182144586620, 0.711075623384, 
                2.161219351339, -1.736085365532, 0.044487219848, 
                0.814387504674, 0.103122826790, -0.831114843672, 
                2.261170992239, -1.094987814727, 0.072489047792])

_x3 = np.array([])

class TestOTPBulk(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.nmol = 5
        self.boxvec = np.array([5,5,5])
        self.rcut = 2.5
        self.system = OTPBulk(self.nmol, self.boxvec, self.rcut)
        pot = self.system.get_potential()
        self.db = self.system.create_database()
        self.m1 = self.db.addMinimum(pot.getEnergy(_x1), _x1)
        self.m2 = self.db.addMinimum(pot.getEnergy(_x2), _x2)
            
    def test_periodic_distance(self):
        x0 = self.system.get_random_configuration()
        shift = np.zeros(self.nmol*6)
        for i in xrange(3*self.nmol):
            shift[i] = self.boxvec[i%3]/2+1
        x1 = x0 + shift
        
        self.assertLess(sqrt(self.system.aatopology.distance_squared(x0,x1)), 
                        np.linalg.norm((self.boxvec/2+1))*len(self.system.aatopology.sites))
 
        shift = np.zeros(self.nmol*6)    
        for i in xrange(3*self.nmol):
            shift[i] = -self.boxvec[i%3]
        x1 = x0 + shift
        
        self.assertLess(self.system.aatopology.distance_squared(x0,x1), 1e-8)     
    
    def test1(self):
        pot = self.system.get_potential()
        self.assertLess(np.linalg.norm(pot.getGradient(self.m1.coords)), .1)
        self.assertLess(np.linalg.norm(pot.getGradient(self.m2.coords)), .1)
    
    def test_basinhopping(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.setPrinting(ostream=None)
        bh.run(5)
        self.assertGreaterEqual(db.number_of_minima(), 1)

    # Not totally convinced the pair of minima I've given will ever be connected by this
    # test. I think this one needs to be reviewed.
#     def test_double_ended_connect(self):
#         connect = self.system.get_double_ended_connect(self.m1, self.m2, self.db)
#         connect.connect()
#         self.assertTrue(connect.success())
             

if __name__ == "__main__":
    unittest.main()

