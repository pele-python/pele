import unittest

import numpy as np

from pele.landscape import DoubleEndedConnect
from pele.transition_states.tests.test_NEB import _x1, _x2


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
        db.addMinimum(pot.getEnergy(_x1), _x1)
        db.addMinimum(pot.getEnergy(_x2), _x2)
        m1, m2 = db.minima()
        
        connect = DoubleEndedConnect(m1, m2, pot, mindist, db, verbosity=10)
        connect.connect()
        self.assertTrue(connect.success())
        
        path = connect.returnPath()

if __name__ == "__main__":
    unittest.main()

