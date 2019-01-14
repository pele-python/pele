from __future__ import print_function
import unittest


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

_x3 = np.array([-0.264936050189,  0.44177611652,  -2.40068725812,
                -0.395208042106, -0.530124165016, -2.00197971185,
                -0.383451054068, -0.0306274737548, 1.29237018543,
                -0.336122426052,  0.0351631353222, 0.179277207348,
                -0.270922039141,  0.133577696802, -0.954901208908,
                 1.02007635412,  -2.4460283088,   -1.63025739996,
                 -1.95972769035,  0.717896436644,  1.04062511399,
                 -1.55924037321, -2.53690688495,   0.973080106314,
                 1.97604673517,  -1.87955962258,   0.658859696552,
                 -0.630973896475, 0.275513983687,  1.97234452197])

_x4 = np.array([-0.277756113289, 0.497817364111, -2.06835022969,
                -0.444788034453, -0.547430061919, -2.12768687134,
                -0.46324808614, 0.313065685544, 1.27172113057,
                -0.215178117231, 0.144054644281, -0.0264648652671,
                -0.249669260443, -0.357742322129, -0.935139950397,
                1.10148972447, -2.27286223172, -2.11965707809,
                -1.59911473012, 0.491623290371, 1.09852365416,
                -1.1150448545, -2.18422190695, 2.34032994013,
                2.02035026408, -2.4866568981, 1.06785759097,
                -0.271404765585, -0.191103624186, 2.00769383374])

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
        self.m3 = self.db.addMinimum(pot.getEnergy(_x3), _x3)
        self.m4 = self.db.addMinimum(pot.getEnergy(_x4), _x4)                
            
    def test_energy(self):
        self.assertAlmostEqual(self.m3.energy, -34.6758148083, 5) 
        self.assertAlmostEqual(self.m4.energy, -33.7606288716, 5)               
            
    def test_periodic_distance(self):
        x0 = self.system.get_random_configuration()
        shift = np.zeros(self.nmol*6)
        for i in range(3*self.nmol):
            shift[i] = self.boxvec[i%3]/2+1
        x1 = x0 + shift
        
        self.assertLess(sqrt(self.system.aatopology.distance_squared(x0,x1)), 
                        np.linalg.norm((self.boxvec//2+1))*len(self.system.aatopology.sites))
 
        shift = np.zeros(self.nmol*6)    
        for i in range(3*self.nmol):
            shift[i] = -self.boxvec[i%3]
        x1 = x0 + shift
        
        self.assertLess(self.system.aatopology.distance_squared(x0,x1), 1e-8)     
    
    def test1(self):
        pot = self.system.get_potential()
        self.assertLess(np.linalg.norm(pot.getGradient(self.m1.coords)), .1)
        self.assertLess(np.linalg.norm(pot.getGradient(self.m2.coords)), .1)
        
    def test_random_configuration(self):
        n_tests = 100
        fail_count = 0
        x = []
        for i in range(n_tests):
            x.append(self.system.get_random_configuration())
            for j in range(i):
                if(i == j):
                    continue
                if (np.linalg.norm(x[i]-x[j]) < 1e-10):
                    fail_count += 1
                    print("Failing configurations:")
                    print(x[i])
                    print(x[j])
                    print("Difference")
                    print(x[i]-x[j])
                    print("Norm")
                    print(np.linalg.norm(x[i]-x[j]))
        if fail_count > 0:
            print("Failed {} times".format(fail_count))
        self.assertEqual(fail_count, 0)
    
    def test_distance_measure(self):
        n_tests = 100
        fail_count = 0
        for i in range(n_tests):
            x1 = self.system.get_random_configuration()
            x2 = self.system.get_random_configuration()
            dist1 = sqrt(self.system.aatopology.distance_squared(x1, x2))
            x1at = self.system.aatopology.to_atomistic(x1)
            x2at = self.system.aatopology.to_atomistic(x2)
            dist2 = np.linalg.norm(x1at - x2at)
            if(dist1-dist2>1e-13):
                fail_count+=1
                print("Failed.")
                print("distance measure:", dist1)
                print("atomistic cartesian distance:", dist2)
        self.assertEqual(fail_count, 0)
            
        
    
    def test_basinhopping(self):
        db = self.system.create_database()
        bh = self.system.get_basinhopping(db)
        bh.setPrinting(ostream=None)
        bh.run(5)
        self.assertGreaterEqual(db.number_of_minima(), 1)

    def test_double_ended_connect(self):
        connect = self.system.get_double_ended_connect(self.m3, self.m4, self.db)
        connect.connect()
        self.assertTrue(connect.success())
             

if __name__ == "__main__":
    unittest.main()

