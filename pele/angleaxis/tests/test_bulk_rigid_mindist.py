from __future__ import print_function
import unittest
import copy
import numpy as np
from math import pi, cos, sin

from pele.angleaxis.rigidbody import RBTopologyBulk, RigidFragmentBulk
from pele.angleaxis.aaperiodicttransforms import MeasurePeriodicRigid, ExactMatchRigidPeriodic, TransformPeriodicRigid
from pele.angleaxis.bulk_rigid_mindist import MinDistBulkRigid

class TestExactMatchPeriodicRigid(unittest.TestCase):
    def setUp(self):
        self.nrigid = 1
        self.boxl = np.array([5,5,5])
      
#         self.topology = RBTopologyBulk(self.boxl)
#         self.topology.finalize_setup()
#         sites = []
#         for i in range(self.nrigid):
#             sites.append(self.make_molecule())
#         self.topology.add_sites(sites)
#         self.draw_bonds = []
#         for i in xrange(self.nrigid):
# #             self.draw_bonds.append((2*i,2*i+1))
#             self.draw_bonds.append((3*i, 3*i+1))
#             self.draw_bonds.append((3*i, 3*i+2))        
#         #self.permlist = self.get_permlist()     
        
        import pele.angleaxis._otp_bulk as OTP        
        self.system = OTP.OTPBulk(self.nrigid, self.boxl, 2.5)
        self.topology = self.system.aatopology
        self.mindist = self.system.get_mindist()       

        self.transform = TransformPeriodicRigid()         
        self.measure = MeasurePeriodicRigid(self.topology, self.transform)
        self.exact_match = ExactMatchRigidPeriodic(self.measure, accuracy=1e-5)
#         self.mindist = MinPermDistBulk(self.boxl, self.measure)      
                
        self.x1 = self.get_random_configuration()
        self.x2diff = self.get_random_configuration()
        self.x2same = self.x1.copy()
        self.x2trans = self.x1.copy()
        trans = np.random.random(3)*self.boxl
        self.transform.translate(self.x2trans, trans)
        #dist = self.measure.get_dist(self.x1, self.x2diff)

                 


    def make_molecule(self):
        """this constructs a single test molecule"""
        molecule = RigidFragmentBulk(self.boxl) 
#         molecule.add_atom("O", np.array([-0.25,0.0,0.0]), 1.)
#         molecule.add_atom("O", np.array([0.25,0.,0.]), 1.)
               
#         Note: need to change setup of self.draw_bonds also.       
               
        molecule.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        molecule.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        molecule.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        molecule.finalize_setup()
        return molecule
    
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
        x = np.zeros([self.nrigid*2,3])
        for i in range(self.nrigid):
            for j in range(3):
                x[i][j] = np.random.uniform(-self.boxl[j]/2., self.boxl[j]/2.)
        for i in range(self.nrigid,2*self.nrigid):
            x[i] = 5.*np.random.random(3)
        #print "x.flatten()", x.flatten()
        return x.flatten()        
    
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2trans))
         
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.topology.boxvec 
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_align_translate(self):
        np.random.seed(4)
        n = self.nrigid
                 
        fail_counter = 0
        for i in range(10000):
#             if (i%100 == 0):
#                 print i
            self.x1 = self.get_random_configuration()
            self.x2diff = self.get_random_configuration()
             
            try:        
                dist, x1, x2 = self.mindist(self.x1, self.x2diff) 
            except RuntimeError:
                pass            
 
            if self.exact_match(self.x2diff,x2) is False:
                fail_counter += 1
                         
        self.assertFalse(fail_counter, "bond lengths were changed %d times" % fail_counter)                
        
    def test_align_match(self):
        np.random.seed(2)        
        fail_counter = 0
        for i in range(1):
#             if (i%100 == 0):
#                 print i
            self.x1 = self.get_random_configuration()
            self.x2trans = self.x1
            translate = np.random.random(3)*self.boxl
            self.transform.translate(self.x2trans, translate)
            try:
                dist, x1, x2 = self.mindist(self.x1, self.x2trans)
            except RuntimeError:
                pass
            if (dist>1e-5):
                fail_counter+=1
                print(dist)
        self.assertFalse(fail_counter, "structure matching failed %d times" % fail_counter)        
        
    def test_align_permutation(self):
        np.random.seed(4)        
        fail_counter = 0
        for i in range(1):
            if (i%100 == 0):
                print(i)
            self.x1 = self.get_random_configuration()
            self.x2diff = self.get_random_configuration()
            try:
                dist2, x1, x2 = self.mindist(self.x1, self.x2diff)    
            except RuntimeError:
                fail_counter+=1
        self.assertFalse(fail_counter, "rotational alignment failed %d times" % fail_counter)                
            
    def test_align_improvement(self, verbose = True):
        np.random.seed(6)
        max_step = 10000        
        fail_counter = 0
        ave = 0
        ave_inc = 0
        
        for i in range(max_step):
            if (i%100 == 0):
                print(i)
            self.x1 = self.get_random_configuration()
            self.x2diff = self.get_random_configuration()
            dist = self.measure.get_dist(self.x1, self.x2diff)
            dist2, x1, x2 = self.mindist(self.x1, self.x2diff)

            if(dist2 > dist):
#             if(i==10):
                fail_counter += 1
                ave_inc += dist2 - dist
                if(verbose):
                    print("x1", x1)
                    print("old x2", self.x2diff)
                    print("new x2", x2)
#                     print "atomistic x1", self.topology.to_atomistic(x1).flatten()
#                     print "atomistic x2", self.topology.to_atomistic(self.x2diff)
#                     print "new atomistic x2", self.topology.to_atomistic(x2)
                    print("dist", dist)
                    print("dist2", dist2)
                    print("i", i)

                try: 
                    import pele.utils.pymolwrapper as pym
                    pym.start()
 
                    x1 = self.topology.to_atomistic(x1)
                    x2 = self.topology.to_atomistic(x2)
                    self.x2diff = self.topology.to_atomistic(self.x2diff)
  
                    pym.draw_rigid(x1, "A", 0.1, (1,0,0), self.system.draw_bonds)
                    pym.draw_rigid(self.x2diff, "B", 0.1, (0,1,0), self.system.draw_bonds)
                    pym.draw_rigid(x2, "C", 0.1, (0,0,1), self.system.draw_bonds) 
                    pym.draw_box(self.boxl, "D", 0.1)  
                    
#                     pym.draw_rigid(x1[:self.nrigid*3], "A", 0.1, (1,0,0))
#                     pym.draw_rigid(self.x2diff[:self.nrigid*3], "B", 0.1, (0,1,0))
#                     pym.draw_rigid(x2[:self.nrigid*3], "C", 0.1, (0,0,1)) 
#                     pym.draw_box(self.boxl, "D", 0.1)                     
                                                                                   
                except:
                    print("Could not draw using pymol, skipping this step")
                    
            else:
                ave += dist - dist2
                
        ave = ave/max_step
        if (fail_counter>0): ave_inc = ave_inc / fail_counter
        print("average decrease in distance", ave)
        print("average increase in distance", ave_inc)

        self.assertFalse(fail_counter, "alignment failed %d times" % fail_counter)


    
if __name__ == '__main__':
    unittest.main()
