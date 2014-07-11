import unittest

import numpy as np
from numpy import cos, sin, pi

from pele.angleaxis import RBTopology, RigidFragment, RBPotentialWrapper
from pele.potentials import LJ

class TestOTP(unittest.TestCase):
    
    def make_otp(self):
        """this constructs a single OTP molecule"""
        otp = RigidFragment()
        otp.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        otp.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        otp.finalize_setup()
        print "otp"
        print otp.atom_positions
        return otp

    
    def setUp(self):
        nrigid = 3
        self.topology = RBTopology()
        self.topology.add_sites([self.make_otp() for i in xrange(nrigid)])
        
        cartesian_potential = LJ()
        self.pot = RBPotentialWrapper(self.topology, cartesian_potential)
        
        self.x0 = [2.550757898788, 2.591553038507, 3.696836364193, 
                2.623281513163, 3.415794212648, 3.310786279789, 
                1.791383852327, 2.264321752809, 4.306217333671, 
                0.761945654023, -0.805817782109, 1.166981882601, 
                0.442065301864, -2.747066418223, -1.784325262714, 
                -1.520905562598, 0.403670860200, -0.729768985400, ]
        self.x0 = np.array(self.x0)
        self.e0 = -17.3387670023
        assert nrigid * 6 == self.x0.size
        
        self.x0atomistic = [ 3.064051819556, 2.474533745459, 3.646107658946,
                            2.412011983074, 2.941152759499, 4.243695098053, 
                            2.176209893734, 2.358972610563, 3.200706335581, 
                            2.786627589565, 3.211876105193, 2.850924310983, 
                            1.962626909252, 3.436918873216, 3.370903763850,
                            3.120590040673, 3.598587659535, 3.710530764535, 
                            1.697360211099, 2.317229950712, 4.823998989452, 
                            2.283487958310, 1.840698306602, 4.168734267290, 
                            1.393303387573, 2.635037001113, 3.925918744272, ]
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
        x0 = np.array(range(self.nrigid * 6), dtype=float);
        x2 = x0.reshape([-1,3])
        for p in x2[self.nrigid:,:]:
            p /= np.linalg.norm(p);
        print x0
        print "range to atomistic"
        print x0
        atomistic = self.topology.to_atomistic(x0).flatten()
        print atomistic
        print atomistic.size
        print atomistic[14]
        print atomistic[23]



if __name__ == "__main__":
    unittest.main()

