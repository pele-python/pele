from __future__ import print_function
import unittest
import copy
import numpy as np

from pele.mindist.periodic_exact_match import MeasurePeriodic, TransformPeriodic, ExactMatchPeriodic
from pele.mindist.periodic_mindist import MinDistBulk
from pele.systems.morse_bulk import MorseBulk

class TestExactMatchPeriodic(unittest.TestCase):
    def setUp(self):
        self.nrigid = 5
        self.boxl = np.array([5,5,5])       
        self.system = MorseBulk(self.nrigid, self.boxl)

#         self.mindist = self.system.get_mindist()       
        self.transform = TransformPeriodic()         
        self.measure = MeasurePeriodic(self.boxl, permlist=[])
        self.exact_match = ExactMatchPeriodic(self.measure, accuracy=1e-5)
        self.mindist = MinDistBulk(self.boxl, self.measure)      
                
        self.x1 = self.system.get_random_configuration()
        self.x2diff = self.system.get_random_configuration()
        self.x2same = self.x1.copy()
        self.x2trans = self.x1.copy()
        trans = np.random.random(3)*self.boxl
        self.transform.translate(self.x2trans, trans)
    
    def randomly_permute(self, x):
        import random
        x = x.reshape(-1,3)
        xnew = x.copy()
        for atomlist in self.permlist:
            permutation = copy.copy(atomlist)
            random.shuffle(permutation)
            xnew[atomlist,:] = x[permutation,:]
        return xnew.flatten()
      
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2trans))
         
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.boxl
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
        
    # This test is pretty superfluous for the atomic-system case. But it doesn't take long and probably
    # doesn't hurt.    
    def test_align_translate(self):
        np.random.seed(4)
        n = self.nrigid
                 
        fail_counter = 0
        for i in range(10000):
#             if (i%100 == 0):
#                 print i
            self.x1 = self.system.get_random_configuration()
            self.x2diff = self.system.get_random_configuration()
             
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
            self.x1 = self.system.get_random_configuration()
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
            self.x1 = self.system.get_random_configuration()
            self.x2diff = self.system.get_random_configuration()
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
            self.x1 = self.system.get_random_configuration()
            self.x2diff = self.system.get_random_configuration()
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

#                 try: 
#                     import pele.utils.pymolwrapper as pym
#                     pym.start()
#  
#                     x1 = self.topology.to_atomistic(x1)
#                     x2 = self.topology.to_atomistic(x2)
#                     self.x2diff = self.topology.to_atomistic(self.x2diff)
#   
#                     pym.draw_rigid(x1, "A", 0.1, (1,0,0), self.system.draw_bonds)
#                     pym.draw_rigid(self.x2diff, "B", 0.1, (0,1,0), self.system.draw_bonds)
#                     pym.draw_rigid(x2, "C", 0.1, (0,0,1), self.system.draw_bonds) 
#                     pym.draw_box(self.boxl, "D", 0.1)  
#                     
# #                     pym.draw_rigid(x1[:self.nrigid*3], "A", 0.1, (1,0,0))
# #                     pym.draw_rigid(self.x2diff[:self.nrigid*3], "B", 0.1, (0,1,0))
# #                     pym.draw_rigid(x2[:self.nrigid*3], "C", 0.1, (0,0,1)) 
# #                     pym.draw_box(self.boxl, "D", 0.1)                     
#                                                                                    
#                 except:
#                     print "Could not draw using pymol, skipping this step"
                    
            else:
                ave += dist - dist2
                
        ave = ave/max_step
        if (fail_counter>0): ave_inc = ave_inc / fail_counter
        print("average decrease in distance", ave)
        print("average increase in distance", ave_inc)

        self.assertFalse(fail_counter, "alignment failed %d times" % fail_counter)


    
if __name__ == '__main__':
    unittest.main()
