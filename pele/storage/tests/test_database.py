import unittest
import os

from pele.storage import Database

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

    def test_add_minimum(self):
        # add a duplicate minimum and ensure the db doesn't change
        self.db.addMinimum(0, [0])
        self.assertEqual(self.db.number_of_minima(), self.nminima)
        

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

    def test_remove_ts(self):
        ts = self.db.transition_states()[0]
        self.db.remove_transition_state(ts)
        self.assertEqual(self.db.number_of_transition_states(), self.nts-1)
        self.assertNotIn(ts, self.db.transition_states())
        
        # m should have 2 minima.  both of those should be gone
        self.assertEqual(self.db.number_of_minima(), self.nminima)

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
    
    def test_minimum_adder_Ecut(self):
        ma = self.db.minimum_adder(Ecut=0)
        n0 = self.db.number_of_minima()
        ma(101., [101.])
        self.assertEqual(n0, self.db.number_of_minima())
        ma(-101., [-101.])
        self.assertEqual(n0+1, self.db.number_of_minima())

    def test_minimum_adder_commit_interval(self):
        ma = self.db.minimum_adder(commit_interval=2)
        # replace db.session.commit with a wrapper that keeps track of how
        # many times it's been called
        self.real_commit = self.db.session.commit
        self.count = 0
        def commit_wrapper():
            self.count += 1
            self.real_commit()
        self.db.session.commit = commit_wrapper
        # commit should be called for the first minimum
        ma(101., [101.])
        self.assertEqual(self.count, 1)
        self.assertEqual(self.nminima+1, self.db.number_of_minima())
        # commit should not be called for the second minimum
        ma(102., [102.])
        self.assertEqual(self.count, 1)
        self.assertEqual(self.nminima+2, self.db.number_of_minima())
        # yes for the third, no for the 4th
        ma(103., [103.])
        self.assertEqual(self.count, 2)
        ma(104., [104.])
        self.assertEqual(self.count, 2)
        # commit should be called when minimum adder is deleted
        del ma
        self.assertEqual(self.count, 3)
        self.assertEqual(self.nminima+4, self.db.number_of_minima())
        
    
    def test_merge_minima(self):
        m1 = self.db.minima()[0]
        m2 = self.db.minima()[1]
        self.db.mergeMinima(m1, m2)
        self.assertEqual(len(self.db.minima()), self.nminima-1)
        # transition states shouldn't be deleted
        self.assertEqual(len(self.db.transition_states()), self.nts)
    
    def test_number_of_minima(self):
        self.assertEqual(self.nminima, self.db.number_of_minima())
    
    def test_number_of_transition_states(self):
        self.assertEqual(self.nts, self.db.number_of_transition_states())
    
    def test_highest_energy_minimum(self):
        m1 = self.db._highest_energy_minimum()
        m2 = self.db.minima()[-1]
        self.assertEqual(m1, m2)
    
    def test_maximum_number_of_minima(self):
        m = self.db.addMinimum(-1., [-1.], max_n_minima=self.nminima)
        self.assertEqual(self.nminima, self.db.number_of_minima())
        self.assertIn(m, self.db.minima())

    def test_maximum_number_of_minima_largestE(self):
        e = float(self.nminima + 1)
        m = self.db.addMinimum(e, [e], max_n_minima=self.nminima)
        self.assertEqual(self.nminima, self.db.number_of_minima())
        self.assertIsNone(m)
        
        #ensure the highest energy minimum is still in the database
        mmax = self.db._highest_energy_minimum()
        self.assertEqual(mmax.energy, float(self.nminima-1))
        
    def test_maximum_number_of_minima_minima_adder(self):
        ma = self.db.minimum_adder(max_n_minima=self.nminima)
        m = ma(-1., [-1.])
        self.assertEqual(self.nminima, self.db.number_of_minima())
        self.assertIn(m, self.db.minima())

    def test_getTSfromID(self):
        ts = self.db.transition_states()[0]
        ts1 = self.db.getTransitionStateFromID(ts._id)
        self.assertEqual(ts, ts1)

    def test_property(self):
        # add some system properties and ensure they were added correctly
        self.db.add_property("natoms", 10)
        p = self.db.get_property("natoms")
        self.assertEqual(p.value(), 10)
        self.db.add_property("eps", 2.1)
        p = self.db.get_property("eps")
        self.assertEqual(p.value(), 2.1)
        self.db.add_property("author", "Jake")
        p = self.db.get_property("author")
        self.assertEqual(p.value(), "Jake")
        self.db.add_property("data", [1, 2])
        p = self.db.get_property("data")
        self.assertEqual(p.value(), [1, 2])
        
        # assert that the not set values are None
        p = self.db.get_property("natoms")
        self.assertIsNone(p.string_value)
        self.assertIsNone(p.float_value)
        self.assertIsNone(p.pickle_value)

        
        p = self.db.get_property("noprop")
        self.assertIsNone(p)

        props = self.db.properties(as_dict=True)
        self.assertIsInstance(props, dict)
        self.assertDictContainsSubset(dict(natoms=10), props)
        self.assertEqual(len(props.items()), 4)
        
        props = self.db.properties(as_dict=False)
        self.assertIn(("natoms", 10), [p.item() for p in props])
        self.assertEqual(len(props), 4)

    def test_property_dtype(self):
        # add some system properties and ensure they were added correctly
        self.db.add_property("natoms", 10, dtype="int")
        p = self.db.get_property("natoms")
        self.assertEqual(p.value(), 10)
        self.db.add_property("eps", 2.1, dtype="float")
        p = self.db.get_property("eps")
        self.assertEqual(p.value(), 2.1)
        self.db.add_property("author", "Jake", dtype="string")
        p = self.db.get_property("author")
        self.assertEqual(p.value(), "Jake")
        self.db.add_property("data", [1, 2], dtype="pickle")
        p = self.db.get_property("data")
        self.assertEqual(p.value(), [1, 2])

    def test_bad_property(self):
        self.db.add_property("natoms", 10, dtype="string")
        p = self.db.get_property("natoms")
        self.assertEqual(p.value(), "10")

    def test_property_overwrite_false(self):
        self.db.add_property("natoms", 10)
        # overwriting with the same value should not raise error
        self.db.add_property("natoms", 10, overwrite=False) 
        # overwriting with a different value should raise error
        with self.assertRaises(RuntimeError):
            self.db.add_property("natoms", 11, overwrite=False)
    
    
    def test_add_properties(self):
        props = dict(natoms=10, author="jake")
        self.db.add_properties(props)
        for name, value in props.iteritems():
            p = self.db.get_property(name)
            self.assertEqual(p.value(), value)
    
    def test_load_wrong_schema(self):
        current_dir = os.path.dirname(__file__)
        dbname = current_dir + "/lj6_schema1.sqlite"
        with self.assertRaises(IOError):
            db = Database(dbname, createdb=False)
    
    def test_load_right_schema(self):
        current_dir = os.path.dirname(__file__)
        dbname = current_dir + "/lj6_schema2.sqlite"
        db = Database(dbname, createdb=False)
    
    def test_invalid(self):
        m = self.db.minima()[0]
        self.assertFalse(m.invalid)
        m.invalid = True
        self.db.session.commit()
        self.assertTrue(m.invalid)
        
        # now with a ts
        m = self.db.transition_states()[0]
        self.assertFalse(m.invalid)
        m.invalid = True
        self.db.session.commit()
        self.assertTrue(m.invalid)
    
    def test_user_data(self):
        m = self.db.minima()[0]
        m.user_data = dict(key="value")
        self.db.session.commit()
        v = m.user_data["key"]
        
        # now with a transition state
        m = self.db.transition_states()[0]
        m.user_data = dict(key="value")
        self.db.session.commit()
        v = m.user_data["key"]
        
def benchmark_number_of_minima():
    import time, sys
    import numpy as np
    db = Database("test.large.db")
    
    if True:
        istart = np.random.randint(0, sys.maxint)
        for i in xrange(istart,istart+10000):
            e = float(i)
            db.addMinimum(e, [e], commit=False)
        db.session.commit()
    else:
        i=1
    
    t1 = time.clock()
    print db.number_of_minima()
    print "time", t1 - time.clock(); t1 = time.clock()
    print db.number_of_minima()
    print "time", t1 - time.clock(); t1 = time.clock()
    print db.number_of_minima()
    print "time", t1 - time.clock(); t1 = time.clock()
    e = float(i+1)
    db.addMinimum(e, [e], commit=False)
    t1 = time.clock()
    print db.number_of_minima()
    print "time", t1 - time.clock(); t1 = time.clock()

    print len(db.minima())
    print "time", t1 - time.clock(); t1 = time.clock()

if __name__ == "__main__":
#    benchmark_number_of_minima()
    unittest.main()