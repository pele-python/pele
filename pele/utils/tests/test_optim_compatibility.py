from __future__ import print_function
import unittest
import os
import hashlib
import tempfile


import numpy as np

from pele.utils.optim_compatibility import OptimDBConverter, WritePathsampleDB
from pele.storage import Database
from pele.potentials.tests import _base_test

class TestOptimCompatibility(unittest.TestCase):
    def test1(self):
        current_dir = os.path.dirname(__file__)
        db = Database()
        converter = OptimDBConverter(db, 
                                     mindata=os.path.join(current_dir, "min.data"),
                                     tsdata=os.path.join(current_dir, "ts.data"),
                                     pointsmin=os.path.join(current_dir, "points.min"),
                                     pointsts=os.path.join(current_dir, "points.ts"),
                                     endianness="<")
        converter.convert()
        self.assertEqual(db.number_of_minima(), 2)
        self.assertEqual(db.number_of_transition_states(), 1)
        
        ts = db.transition_states()[0]
        self.assertAlmostEqual(ts.coords[0], 1.548324, 5)
        self.assertAlmostEqual(ts.coords[2], 0.18178001, 5)
        self.assertAlmostEqual(ts.coords[-1], -0.50953229, 5)
        self.assertEqual((180,), ts.coords.shape)

        m = db.minima()[0]
        print(repr(m.coords))
        self.assertAlmostEqual(m.coords[0], 1.53700142, 5)
        self.assertAlmostEqual(m.coords[2], 0.87783657, 5)
        self.assertAlmostEqual(m.coords[-1], -0.50953229, 5)
        self.assertEqual((180,), m.coords.shape)
    
    def test_nocoords(self):
        current_dir = os.path.dirname(__file__)
        db = Database()
        converter = OptimDBConverter(db, 
                                     mindata=os.path.join(current_dir, "min.data"),
                                     tsdata=os.path.join(current_dir, "ts.data"),
                                     pointsmin="poo",
                                     pointsts="poo",
                                     assert_coords=False,
                                     endianness="<")
        converter.convert()

    def test_nocoords_Fails(self):
        current_dir = os.path.dirname(__file__)
        db = Database()
        converter = OptimDBConverter(db, 
                                     mindata=os.path.join(current_dir, "min.data"),
                                     tsdata=os.path.join(current_dir, "ts.data"),
                                     pointsmin="poo",
                                     pointsts="poo",
                                     assert_coords=True,
                                     endianness="<")
        with self.assertRaises(IOError):
            converter.convert()

def sha1_of_file(filename):
    h = hashlib.sha1()
    with open(filename, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


class TestWritePathsampleDB(unittest.TestCase):
    
    def test1(self):
        current_dir = os.path.dirname(__file__)
        db = Database()
        converter = OptimDBConverter(db, 
                                     mindata=os.path.join(current_dir, "min.data"),
                                     tsdata=os.path.join(current_dir, "ts.data"),
                                     pointsmin=os.path.join(current_dir, "points.min"),
                                     pointsts=os.path.join(current_dir, "points.ts"),
                                     endianness="<")
        
        converter.convert()
        
        mdata = tempfile.NamedTemporaryFile(delete=True)
        tsdata = tempfile.NamedTemporaryFile(delete=True)
        pm = tempfile.NamedTemporaryFile(delete=True)
        pts = tempfile.NamedTemporaryFile(delete=True)
        print(mdata.name, tsdata.name)
        writer = WritePathsampleDB(db,
                                   mindata=mdata.name,
                                   tsdata=tsdata.name,
                                   pointsmin=pm.name,
                                   pointsts=pts.name,
                                   endianness="<")
        
        writer.write_db()
        
        d1 = np.genfromtxt(os.path.join(current_dir, "min.data"))[:,:3]
        d2 = np.genfromtxt(mdata.name)[:,:3]
        _base_test.assert_arrays_almost_equal(self, d1, d2)
        
        d1 = np.genfromtxt(os.path.join(current_dir, "ts.data")).reshape(-1,8)[:,:5]
        d2 = np.genfromtxt(tsdata.name).reshape(-1,8)[:,:5]
        print(d1, d2)
        _base_test.assert_arrays_almost_equal(self, d1, d2)
        
        self.assertEqual(sha1_of_file(os.path.join(current_dir, "points.min")),
                         sha1_of_file(pm.name))
        self.assertEqual(sha1_of_file(os.path.join(current_dir, "points.ts")),
                         sha1_of_file(pts.name))
        
    def test2(self):
        from pele.utils.tests.test_disconnectivity_graph import create_random_database
        
        db = create_random_database(nmin=20, nts=15)
        
        delete=False
        mdata = tempfile.NamedTemporaryFile(delete=delete, suffix=".min.data")
        tsdata = tempfile.NamedTemporaryFile(delete=delete, suffix=".ts.data")
        pm = tempfile.NamedTemporaryFile(delete=delete)
        pts = tempfile.NamedTemporaryFile(delete=delete)
        print(mdata.name, tsdata.name)
        writer = WritePathsampleDB(db,
                                   mindata=mdata.name,
                                   tsdata=tsdata.name,
                                   pointsmin=pm.name,
                                   pointsts=pts.name,
                                   endianness="<")
        
        writer.write_db()
        
        newdb = Database()
        reader = OptimDBConverter(newdb,
                                  mindata=mdata.name,
                                  tsdata=tsdata.name,
                                  pointsmin=pm.name,
                                  pointsts=pts.name,
                                  endianness="<")
        reader.convert()
        
        for m1, m2 in zip(db.minima(), newdb.minima()):
            self.assertAlmostEqual(m1.energy, m2.energy)
            _base_test.assert_arrays_almost_equal(self, m1.coords, m2.coords)
        
        for ts1, ts2 in zip(db.transition_states(), newdb.transition_states()):
            self.assertAlmostEqual(ts1.energy, ts2.energy)
            _base_test.assert_arrays_almost_equal(self, ts1.coords, ts2.coords)
            self.assertAlmostEqual(ts1.minimum1.energy, ts2.minimum1.energy)
            self.assertAlmostEqual(ts1.minimum2.energy, ts2.minimum2.energy)
        
        
                                   

if __name__ == "__main__":
    unittest.main()

