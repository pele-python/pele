import unittest

from pele.systems import LJCluster
from pele.takestep import RandomDisplacement
from pele.takestep import GroupSteps, BlockMoves, Reseeding


class RDWrap(RandomDisplacement):
    count = 0
    def takeStep(self, coords, **kwargs):
        self.count += 1 
        RandomDisplacement.takeStep(self, coords, **kwargs)
        
class TestGroupTakestep(unittest.TestCase):
    def test1(self):
        natoms = 13
        system = LJCluster(natoms)
        ts1 = RDWrap()
        ts2 = RDWrap()
        
        ts = GroupSteps([ts1, ts2])
        bh = system.get_basinhopping(takestep=ts)
        bh.run(2)
        self.assertEqual(ts1.count, 2) 
        self.assertEqual(ts2.count, 2) 

class TestBlockTakestep(unittest.TestCase):
    def test1(self):
        natoms = 13
        system = LJCluster(natoms)
        ts1 = RDWrap()
        ts2 = RDWrap()
        
        ts = BlockMoves()
        n1 = 3
        n2 = 5
        ts.addBlock(n1, ts1)
        ts.addBlock(n2, ts2)
        
        bh = system.get_basinhopping(takestep=ts)
        bh.run(n1)
        self.assertEqual(ts1.count, n1) 
        self.assertEqual(ts2.count, 0) 
        
        bh.run(n2)
        self.assertEqual(ts1.count, n1) 
        self.assertEqual(ts2.count, n2)
        
        bh.run(n1 + n2)
        self.assertEqual(ts1.count, 2*n1) 
        self.assertEqual(ts2.count, 2*n2)
        
        bh.run(5*(n1 + n2) - 1)
        self.assertEqual(ts1.count, 7*n1) 
        self.assertEqual(ts2.count, 7*n2 - 1)
    
class TestReseeding(unittest.TestCase):
    def test1(self):
        # lj cluster with 4 atoms has only one minimum
        natoms = 4
        system = LJCluster(natoms)
        
        ts1 = RDWrap()
        reseed = RDWrap()
        maxnoimprove = 20
        ts = Reseeding(ts1, reseed, maxnoimprove=maxnoimprove)
        
        bh = system.get_basinhopping(takestep=ts)
        bh.run(1) # find the global minimum
        bh.run(maxnoimprove - 1)
        
        self.assertEqual(reseed.count, 0)
        self.assertEqual(ts1.count, maxnoimprove)
        
        bh.run(1)
        self.assertEqual(reseed.count, 1)
        self.assertEqual(ts1.count, maxnoimprove)
        
        bh.run(maxnoimprove-1)
        self.assertEqual(reseed.count, 1)
        self.assertEqual(ts1.count, 2*maxnoimprove-1)
        
        bh.run(1)
        self.assertEqual(reseed.count, 2)
        self.assertEqual(ts1.count, 2*maxnoimprove-1)
        
        

if __name__ == "__main__":
    unittest.main()

