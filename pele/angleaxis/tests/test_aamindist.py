import numpy as np
from copy import deepcopy
from pele.angleaxis.molecules import create_water
from pele.angleaxis import RBTopology
import pele.angleaxis.aamindist as am
#import gmin_ as GMIN
from pele.potentials import GMINPotential
import unittest

class TestAAMindist(unittest.TestCase):
    def setUp(self):
        #GMIN.initialize()
#        self.pot = GMINPotential(GMIN)
        self.nrigid = 10
        self.water = create_water()
        self.topology = RBTopology()
        self.topology.add_sites([deepcopy(self.water) for i in xrange(self.nrigid)])    

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
        for i in xrange(100):
            coords1 = np.random.random(6*self.nrigid)*4
            coords2 = np.random.random(6*self.nrigid)*4
            
            coords1[:3*self.nrigid]=0
            coords2[:3*self.nrigid]=0
            
            measure1 = am.MeasureAngleAxisCluster(self.topology)
            measure2 = am.MeasureRigidBodyCluster(self.topology)
            
            #print 
            self.assertAlmostEqual(measure1.get_dist(coords1, coords2),
                                   measure2.get_dist(coords1, coords2))

    
if __name__ == '__main__':
    unittest.main()
