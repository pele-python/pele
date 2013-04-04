from pygmin.storage import Database
import unittest

class TestDB(unittest.TestCase):
    def setUp(self):
        self.db = Database()
        self.nminima = 10
        for i in range(self.nminima):
            e = float(i)
            self.db.addMinimum(e, [e])
        
        
        self.nts = 3
        self.db.addTransitionState(0., [0.], self.db.minima()[0], self.db.minima()[1], eigenval=0., eigenvec=[0.])
        self.db.addTransitionState(0., [0.], self.db.minima()[1], self.db.minima()[2], eigenval=0., eigenvec=[0.])
        self.db.addTransitionState(0., [0.], self.db.minima()[0], self.db.minima()[2], eigenval=0., eigenvec=[0.])

    def test_size(self):
        self.assertEqual(len(self.db.minima()), self.nminima)
        
    def test_energy(self):
        m = self.db.minima()[0]
        self.assertEqual(m.energy, 0.)

    def test_coords(self):
        m = self.db.minima()[0]
        self.assertEqual(m.coords, [0.])

    def test_sizets(self):
        self.assertEqual(len(self.db.transition_states()), self.nts)
    def test_energyts(self):
        ts = self.db.transition_states()[0]
        self.assertEqual(ts.energy, 0.)

    def test_coordsts(self):
        ts = self.db.transition_states()[0]
        self.assertEqual(ts.coords, [0.])
    
    def test_remove_minimum(self):
        m = self.db.minima()[0]
        self.db.removeMinimum(m)
        self.assertEqual(len(self.db.minima()), self.nminima-1)
        self.assertNotIn(m, self.db.minima())
        
        # m should have 2 minima.  both of those should be gone
        self.assertEqual(len(self.db.transition_states()), self.nts-2)

    def test_getTransitionState(self):
        m1 = self.db.minima()[0]
        m2 = self.db.minima()[1]
        m3 = self.db.minima()[-1]
        self.assertIsNotNone(self.db.getTransitionState(m1, m2))
        self.assertIsNone(self.db.getTransitionState(m1, m3))
    
    def test_getMinimum(self):
        m = self.db.minima()[0]
        self.assertEqual(m, self.db.getMinimum(m._id))
        
    def test_minimum_adder(self):
        ma = self.db.minimum_adder()
        ma(101., [101.])
        self.assertEqual(len(self.db.minima()), self.nminima+1)
    
    def test_merge_minima(self):
        m1 = self.db.minima()[0]
        m2 = self.db.minima()[1]
        self.db.mergeMinima(m1, m2)
        self.assertEqual(len(self.db.minima()), self.nminima-1)
        # transition states shouldn't be deleted
        self.assertEqual(len(self.db.transition_states()), self.nts)
        



if __name__ == "__main__":
    unittest.main()