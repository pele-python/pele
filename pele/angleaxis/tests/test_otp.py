import unittest
from itertools import izip

import numpy as np
from numpy import cos, sin, pi

from pele.angleaxis import RBTopology, RigidFragment, RBPotentialWrapper
from pele.potentials import LJ
from pele.angleaxis._otp_cluster import OTPCluster
from pele.thermodynamics import get_thermodynamic_information
from pele.utils import rotations
from pele.angleaxis._aa_utils import _rot_mat_derivative, _sitedist_grad, _sitedist
from pele.angleaxis.aamindist import MeasureRigidBodyCluster


_x03 = np.array([2.550757898788, 2.591553038507, 3.696836364193, 
                2.623281513163, 3.415794212648, 3.310786279789, 
                1.791383852327, 2.264321752809, 4.306217333671, 
                0.761945654023, -0.805817782109, 1.166981882601, 
                0.442065301864, -2.747066418223, -1.784325262714, 
                -1.520905562598, 0.403670860200, -0.729768985400])
_x03_atomistic = np.array([3.064051819556, 2.474533745459, 3.646107658946,
                            2.412011983074, 2.941152759499, 4.243695098053, 
                            2.176209893734, 2.358972610563, 3.200706335581, 
                            2.786627589565, 3.211876105193, 2.850924310983, 
                            1.962626909252, 3.436918873216, 3.370903763850,
                            3.120590040673, 3.598587659535, 3.710530764535, 
                            1.697360211099, 2.317229950712, 4.823998989452, 
                            2.283487958310, 1.840698306602, 4.168734267290, 
                            1.393303387573, 2.635037001113, 3.925918744272
                           ])

class TestOTPExplicit(unittest.TestCase):
    
    def make_otp(self):
        """this constructs a single OTP molecule"""
        otp = RigidFragment()
        otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        otp.finalize_setup()
        return otp

    
    def setUp(self):
        nrigid = 3
        self.topology = RBTopology()
        self.topology.add_sites([self.make_otp() for i in xrange(nrigid)])
        self.topology.finalize_setup()
        
        cartesian_potential = LJ()
        self.pot = RBPotentialWrapper(self.topology, cartesian_potential)
        
        self.x0 = _x03
        self.x0 = np.array(self.x0)
        self.e0 = -17.3387670023
        assert nrigid * 6 == self.x0.size
        
        self.x0atomistic = _x03_atomistic
        self.nrigid = nrigid
    
    def test_energy(self):
        e = self.pot.getEnergy(self.x0)
        self.assertAlmostEqual(e, self.e0, delta=1e-4)

    def test_energy_gradient(self):
        e = self.pot.getEnergy(self.x0)
        gnum = self.pot.NumericalDerivative(self.x0)
         
        e2, g = self.pot.getEnergyGradient(self.x0)
        self.assertAlmostEqual(e, e2, delta=1e-4)
         
        for i in xrange(g.size):
            self.assertAlmostEqual(g[i], gnum[i], 2)
    
    def test_to_atomistic(self):
        xatom = self.topology.to_atomistic(self.x0).flatten()
        for i in xrange(xatom.size):
            self.assertAlmostEqual(xatom[i], self.x0atomistic[i], 2)
    
    def test_site_to_atomistic(self):
        rf = self.make_otp()
        p = np.array([1., 2, 3])
        p /= np.linalg.norm(p)
        com = np.array([4., 5, 6])
        print "otp to atomistic"
        print rf.to_atomistic(com, p)
        

        print "otp transform grad"
        g = np.array(range(9), dtype=float).reshape([-1,3])
        print g.reshape(-1)
        
        print rf.transform_grad(p, g)
    
    def test_to_atomistic2(self):
        x0 = np.array(range(self.nrigid * 6), dtype=float)
        x2 = x0.reshape([-1,3])
        for p in x2[self.nrigid:,:]:
            p /= np.linalg.norm(p)
        atomistic = self.topology.to_atomistic(x0).flatten()
        
        from pele.potentials import LJ
        lj = LJ()
        e, g = lj.getEnergyGradient(atomistic.reshape(-1))
        grb = self.topology.transform_gradient(x0, g)
        rbpot = RBPotentialWrapper(self.topology, lj)
        print rbpot.getEnergy(x0)


class TestCppRBPotentialWrapper(TestOTPExplicit):
    def test_pot_wrapper(self):
        from pele.angleaxis import _cpp_aa
        from pele.potentials import LJ
        rbpot_cpp = _cpp_aa.RBPotentialWrapper(self.topology, LJ())
        rbpot = RBPotentialWrapper(self.topology, LJ())
        
        self.assertAlmostEqual(rbpot_cpp.getEnergy(self.x0), 
                               rbpot.getEnergy(self.x0), 4)
        
        e1, grad1 = rbpot_cpp.getEnergyGradient(self.x0)
        e2, grad2 = rbpot.getEnergyGradient(self.x0)
        self.assertAlmostEqual(e1, e2, 4)
        for g1, g2 in zip(grad1, grad2):
            self.assertAlmostEqual(g1, g2, 3) 
#         print "energy cpp"
#         print e1, e2
#         print grad1
#         print grad2
        

_x1 = np.array([ 1.9025655 ,  0.39575842,  2.70994994,  1.12711741,  0.63413933,
                1.99433564,  1.86553644,  1.71434811,  2.22927686,  0.80189315,
                1.19513512,  3.02357997,  1.25845172, -0.06244027,  1.27217385,
               -2.26564485,  0.25537024,  0.66231258, -1.49510664,  0.94428774,
               -0.04120075, -0.87664883, -0.21441754,  2.05796547])
_x2 = np.array([ 2.01932983,  0.32928065,  2.34949584,  1.12261277,  0.84195098,
                2.08827517,  1.42644916,  1.83608794,  2.23147536,  1.12872074,
                0.93206141,  3.28789605,  1.73243138, -0.1199651 ,  1.02925229,
               -1.64603729,  0.30701482,  0.90204992, -1.96259809,  0.06557119,
                0.11010908, -0.37462588, -0.42374544,  1.97728056])
 
class TestOTPCluster(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.nmol = 4
        self.system = OTPCluster(self.nmol)
        pot = self.system.get_potential()
        self.db = self.system.create_database()
        self.m1 = self.db.addMinimum(pot.getEnergy(_x1), _x1)
        self.m2 = self.db.addMinimum(pot.getEnergy(_x2), _x2)
    
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

    def test_double_ended_connect(self):
        connect = self.system.get_double_ended_connect(self.m1, self.m2, self.db)
        connect.connect()
        self.assertTrue(connect.success())
        
        path = connect.returnPath()
    
    def test_thermodynamics(self):
        get_thermodynamic_information(self.system, self.db, nproc=None, recalculate=True)
        self.assertIsNotNone(self.m1.fvib)
        
        mt = self.system.get_metric_tensor(self.m1.coords)
        print "metric tensor"
        print mt
    
class TestRBTopologyOTP(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)
        self.nmol = 3
        self.system = OTPCluster(self.nmol)
#        pot = self.system.get_potential()
#        self.db = self.system.create_database()
#        self.m1 = self.db.addMinimum(pot.getEnergy(_x1), _x1)
#        self.m2 = self.db.addMinimum(pot.getEnergy(_x2), _x2)

        self.x0 = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 
                             0.517892, 0.575435, 0.632979, 
                             0.531891, 0.576215, 0.620539, 
                             0.540562, 0.5766, 0.612637 ])
        
        from pele.angleaxis.aamindist import TransformAngleAxisCluster
        self.topology = self.system.aatopology
        self.transform = TransformAngleAxisCluster(self.topology)
        
        self.p0 = np.array(range(1,4), dtype=float)
        self.p0 /= np.linalg.norm(self.p0)

    
    def test_transform_rotate(self):
        print "\ntest rotate"
        x = self.x0.copy()
        p = np.array(range(1,4), dtype=float)
        p /= np.linalg.norm(p)
        self.transform.rotate(x, rotations.aa2mx(p))
        
        xnewtrue = np.array([ 0.48757698,  0.61588594,  2.09355038,  2.02484605,  4.76822812,
                            4.81289924,  3.56211511,  8.92057031,  7.53224809,  0.71469473,
                            1.23875927,  1.36136748,  0.72426504,  1.24674367,  1.34426835,
                            0.73015833,  1.25159032,  1.33345003])
        for v1, v2 in izip(x, xnewtrue):
            self.assertAlmostEqual(v1, v2, 5)
    
    def test_align_path(self):
        print "\ntest align_path"
        x1 = self.x0.copy()
        x2 = self.x0 + 5
        
        self.topology.align_path([x1, x2])
        
        x2true = np.array([  5.        ,   6.        ,   7.        ,   8.        ,
                             9.        ,  10.        ,  11.        ,  12.        ,
                            13.        ,   1.92786071,   1.94796529,   1.96807021,
                             1.93320298,   1.94869267,   1.96418236,   1.93645608,
                             1.94905155,   1.96164668])
        
        for v1, v2 in izip(x1, self.x0):
            self.assertAlmostEqual(v1, v2, 5)
        for v1, v2 in izip(x2, x2true):
            self.assertAlmostEqual(v1, v2, 5)
    
    def test_cpp_zero_ev(self):
        print "\ntest zeroEV cpp"
        x = self.x0.copy()
        zev = self.topology._zeroEV_python(x)
        czev = self.topology.cpp_topology.get_zero_modes(x)
        self.assertEqual(len(czev), 6)
        for ev, cev in izip(zev, czev):
            for v1, v2 in izip(ev, cev):
                self.assertAlmostEqual(v1, v2, 5)     
    
    def test_site_distance_squared(self):
        print "\ntest site distance squared"
        c0 = np.zeros(3)
        c1 = np.ones(3)
        p0 = self.p0.copy()
        p1 = p0 + 1
        site = self.system.make_otp()
        d2 = site.distance_squared(c0, p0, c1, p1)
        d2p = _sitedist(c1-c0, p0, p1, site.S, site.W, site.cog)
        self.assertAlmostEqual(d2, 10.9548367929, 5)


    def test_distance_squared(self):
        print "\ntest distance squared"
        x1 = self.x0.copy()
        x2 = self.x0 + 1.1
        d2 = self.topology.distance_squared(x1, x2)
        d3 = self.topology._distance_squared_python(x1, x2)
        self.assertAlmostEqual(d2, 38.9401810973, 5)
        self.assertAlmostEqual(d2, d3, 5)
        


    def test_distance_squared_grad(self):
        print "\ntest distance squared grad"
        x1 = self.x0.copy()
        x2 = self.x0 + 1.1
        grad = self.topology.distance_squared_grad(x1, x2)
        g2 = self.topology._distance_squared_grad_python(x1, x2)
        
        gtrue = np.array([-6.6       , -6.6       , -6.6       , -6.6       , -6.6       ,
                       -6.6       , -6.6       , -6.6       , -6.6       , -1.21579025,
                       -0.07013805, -1.2988823 , -1.21331786, -0.06984532, -1.28945301,
                       -1.2116105 , -0.06975828, -1.28362943])
        for v1, v2 in izip(grad, gtrue):
            self.assertAlmostEqual(v1, v2, 5)
        for v1, v2 in izip(grad, g2):
            self.assertAlmostEqual(v1, v2, 5)
    
    def test_measure_align(self):
        print "\ntest measure align"
        x1 = self.x0.copy()
        x2 = self.x0 + 5.1
        x2[-1] = x1[-1] + .1
        x20 = x2.copy()
        measure = MeasureRigidBodyCluster(self.topology)
        measure.align(x1, x2)

if __name__ == "__main__":
    unittest.main()

