import os
import tempfile
import unittest

import numpy as np

from pele.storage.savenlowest import SaveN


class TestSaveN(unittest.TestCase):
    def test_keeps_lowest(self):
        saven = SaveN(nsave=2)
        for energy in (3.0, 1.0, 2.0, 1.0 + 1e-4):
            saven(energy, np.full(3, energy))
        self.assertEqual([m.energy for m in saven.data], [1.0, 2.0])

    def test_save_load(self):
        saven = SaveN(nsave=2)
        saven(1.0, np.arange(3.0))
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, "saven.pickle")
            saven.save(filename)
            loaded = SaveN.load(filename)
        self.assertEqual([m.energy for m in loaded.data], [1.0])
        np.testing.assert_array_equal(loaded.data[0].coords, np.arange(3.0))
        loaded(0.5, np.zeros(3))  # the lock is restored
        self.assertEqual(len(loaded.data), 2)


if __name__ == "__main__":
    unittest.main()
