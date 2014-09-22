import unittest
import os

import numpy as np

from pele.utils.optim_compatibility import OptimDBConverter
from pele.storage import Database

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
        print repr(m.coords)
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

if __name__ == "__main__":
    unittest.main()
