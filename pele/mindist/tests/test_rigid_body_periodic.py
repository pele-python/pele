import unittest
import copy
import numpy as np
#from math import pi, cos, sin

from pele.mindist.periodic_exact_match import MeasurePeriodicRigid, ExactMatchRigidPeriodic, TransformPeriodicRigid
from pele.angleaxis.rigidbody import RBTopologyBulk, RigidFragmentBulk
from pele.mindist.minpermdist_stochastic import MinPermDistBulk


class TestExactMatchPeriodicRigid(unittest.TestCase):
    def setUp(self):
        self.nrigid = 3
        self.boxl = np.array([15,10,5])
 
        #boxl = (float(self.natoms) / rho)**(1./3)
        #boxlengths = np.ones(3) * boxl + np.random.rand(3)*.1
       
        self.topology = RBTopologyBulk(self.boxl)
        sites = []
        for i in range(self.nrigid):
            sites.append(self.make_molecule())
        self.topology.add_sites(sites)
        #self.permlist = self.get_permlist()      
        self.measure = MeasurePeriodicRigid(self.boxl, self.topology)
        self.transform = TransformPeriodicRigid()
        self.exact_match = ExactMatchRigidPeriodic(self.measure, accuracy=1e-5)
        self.mindist = MinPermDistBulk(self.boxl, self.measure)
                
        self.x1 = self.get_random_configuration()
        self.x2diff = self.get_random_configuration()
        self.x2same = self.x1.copy()
        self.x2trans = self.x1.copy()
        trans = np.random.random(3)*self.boxl
        self.transform.translate(self.x2trans, trans)
        #dist = self.measure.get_dist(self.x1, self.x2diff)

                 
#     def get_permlist(self):
#         return [range(self.natoms)]

    def make_molecule(self):
        """this constructs a single test molecule"""
        molecule = RigidFragmentBulk(self.boxl) 
        molecule.add_atom("O", np.array([-0.25,0.0,0.0]), 1.)
        molecule.add_atom("O", np.array([0.25,0.,0.]), 1.)
               
        #molecule.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        #molecule.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        #molecule.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
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
        x = np.zeros([self.nrigid,6])
        for i in range(3):
            x[:,i] = np.random.uniform(-self.boxl[i]/2., self.boxl[i]/2., self.nrigid)
        for i in range(3,6):
            x[:,i] = 5.*np.random.random(self.nrigid)
        return x.flatten()        
    
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2trans))
        
        
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.measure.boxlengths  
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_align_improvement(self, verbose = True):
        np.random.seed(2)        
        fail_counter = 0
        for i in range(750):
            self.x1 = self.get_random_configuration()
            self.x2diff = self.get_random_configuration()
            dist = self.measure.get_dist(self.x1, self.x2diff)
            dist2, x1, x2 = self.mindist(self.x1, self.x2diff)
            if(dist2 > dist):
                fail_counter += 1
                if(verbose):
                    print "x1", x1
                    print "x2", self.x2diff
                    print "new x2", x2
                    print "dist", dist
                    print "dist2", dist2
                    
                    try: 
                        import pele.utils.pymolwrapper as pym
                        pym.start()
                        pym.draw_rigid(x1, "A", 1, (1,0,0))
                        pym.draw_rigid(x2, "B", 1, (0,1,0))                        
                    except:
                        print "Could not draw using pymol, skipping this step"

        self.assertFalse(fail_counter, "alignment failed %d times" % fail_counter)

#class TestExactMatchPeriodicBLJ(TestExactMatchPeriodicLJ):        
#    def get_permlist(self):
#        self.ntypeA = int(self.natoms *.8)
#        return [range(self.ntypeA), range(self.ntypeA, self.natoms)]
    
if __name__ == '__main__':
    unittest.main()