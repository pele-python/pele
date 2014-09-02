import unittest

import numpy as np

from pele.landscape import DoubleEndedConnect

from pele.potentials import LJ

x1 = np.array([ 1.33553771, -0.68037346,  0.34217452,  0.99605849,  0.3991415 ,
        0.22669776, -0.40689511, -0.69241725, -0.25056957,  1.31739168,
       -1.16646925, -0.68608386,  0.76810225,  0.58022255, -0.87292919,
        0.96669726, -0.38737814, -1.43705923,  0.24050881, -1.17391284,
       -1.05241976,  0.61828363, -0.38738569, -0.41286104,  0.46846513,
       -1.35499386,  0.0472071 , -0.08082435,  0.3916979 , -0.1396383 ,
       -0.09897039, -0.09439797, -1.16789673,  0.26987003, -0.38739314,
        0.61133724,  1.64346237, -0.08235405, -0.57515248])
x2 = np.array([-0.27549725,  1.04951451, -0.67460502,  0.75753839, -0.46638045,
       -1.95958793, -0.15298401,  0.57974577,  0.35842636, -0.65567666,
       -0.44390788,  0.34193545,  0.29068953,  0.3927537 , -1.4215933 ,
        0.46821022, -0.37270016,  0.25937497, -1.1225271 ,  0.41522655,
        0.87992888, -0.08949013, -1.10066902, -0.40505256,  0.75078783,
       -0.48710754, -0.84378191, -0.21200374, -0.63090045, -1.43808399,
       -0.83319699,  0.32154552, -1.33903316, -1.11577523,  0.43595228,
       -0.23587694, -0.18249336, -0.0255772 , -0.53982891])

class TestDoubleEndedConnect(unittest.TestCase):
#    def test1(self):
#        from pele.systems import LJCluster
#        natoms = 13
#        system = LJCluster(natoms)
#        
#        db = system.create_database()
#        
#        # get some minima
#        bh = system.get_basinhopping(database=db, outstream=None)
#        bh.run(100)
#        
#        m1 = db.minima()[0]
#        m2 = db.minima()[-1]
#        print repr(m1.coords)
#        print repr(m2.coords)
#        
#        connect = DoubleEndedConnect(m1, m2, system.get_potential(), system.get_mindist(), db)
#        connect.connect()

    
    def test2(self):
        from pele.storage import Database
        from pele.systems import LJCluster
        np.random.seed(0)

        natoms = 13
        system = LJCluster(natoms)
        pot = system.get_potential()
        mindist = system.get_mindist(niter=1)
        
        db = Database()
        db.addMinimum(pot.getEnergy(x1), x1)
        db.addMinimum(pot.getEnergy(x2), x2)
        m1, m2 = db.minima()
        
        connect = DoubleEndedConnect(m1, m2, pot, mindist, db, verbosity=10)
        connect.connect()
        self.assertTrue(connect.success())
        
        path = connect.returnPath()

if __name__ == "__main__":
    unittest.main()
