import unittest
import math
import numpy as np
import pele.exceptions as exc
from pele.accept_tests.spherical_container import SphericalContainer as sphere


class TestSphericalContainer(unittest.TestCase):
    def setUp(self):
        """
        creates a set of coordinates which pass/fail a spherical container test
        and corresponding spherical container tests
        """
        # Coordinates
        self.coords = [np.random.uniform(-1.0, 1.0) / math.sqrt(3) for _ in range(30)]
        self.coords_0_1 = [x * 0.1 for x in self.coords]
        self.coords_10 = [x * 10 for x in self.coords]
        self.coords_20 = [x * 20 for x in self.coords]
        self.coords_100 = [x * 100 for x in self.coords]
        # Containers
        self.sphere_0 = sphere(0)
        self.sphere_0_1 = sphere(0.1)
        self.sphere_1 = sphere(1.0)
        self.sphere_10 = sphere(10.0)
        self.sphere_20 = sphere(20.0)
        self.sphere_100 = sphere(100.0)

        # containers nocenter
        self.sphere_1_nc = sphere(1.0, nocenter=True)
        self.sphere_20_nc = sphere(20.0, nocenter=True)

    def test_accept(self):
        # test that coordinates with radius smaller than the container are accepted
        # Coords radius 0.1
        self.assertTrue(self.sphere_0_1.accept(self.coords_0_1))
        self.assertTrue(self.sphere_1.accept(self.coords_0_1))
        self.assertTrue(self.sphere_10.accept(self.coords_0_1))
        self.assertTrue(self.sphere_20.accept(self.coords_0_1))
        self.assertTrue(self.sphere_100.accept(self.coords_0_1))
        self.assertTrue(self.sphere_1_nc.accept(self.coords_0_1))
        self.assertTrue(self.sphere_20_nc.accept(self.coords_0_1))
        # Coords radius 1
        self.assertTrue(self.sphere_1.accept(self.coords))
        self.assertTrue(self.sphere_10.accept(self.coords))
        self.assertTrue(self.sphere_20.accept(self.coords))
        self.assertTrue(self.sphere_100.accept(self.coords))
        self.assertTrue(self.sphere_1_nc.accept(self.coords))
        self.assertTrue(self.sphere_20_nc.accept(self.coords))
        # Coords radius 10
        self.assertTrue(self.sphere_10.accept(self.coords_10))
        self.assertTrue(self.sphere_20.accept(self.coords_10))
        self.assertTrue(self.sphere_100.accept(self.coords_10))
        self.assertTrue(self.sphere_20_nc.accept(self.coords_10))
        # Coords radius 20
        self.assertTrue(self.sphere_20.accept(self.coords_20))
        self.assertTrue(self.sphere_100.accept(self.coords_20))
        self.assertTrue(self.sphere_20_nc.accept(self.coords_20))
        # Coords radius 100
        self.assertTrue(self.sphere_100.accept(self.coords_100))

    def test_reject(self):
        # test that coordinates with radius larger than the container are rejected
        # Coords radius 1
        self.assertFalse(self.sphere_0_1.accept(self.coords))
        # Coords radius 10
        self.assertFalse(self.sphere_0_1.accept(self.coords_10))
        self.assertFalse(self.sphere_1.accept(self.coords_10))
        self.assertFalse(self.sphere_1_nc.accept(self.coords_10))
        # Coords radius 20
        self.assertFalse(self.sphere_0_1.accept(self.coords_20))
        self.assertFalse(self.sphere_1.accept(self.coords_20))
        self.assertFalse(self.sphere_10.accept(self.coords_20))
        self.assertFalse(self.sphere_1_nc.accept(self.coords_20))
        # Coords radius 100
        self.assertFalse(self.sphere_0_1.accept(self.coords_100))
        self.assertFalse(self.sphere_1.accept(self.coords_100))
        self.assertFalse(self.sphere_10.accept(self.coords_100))
        self.assertFalse(self.sphere_20.accept(self.coords_100))
        self.assertFalse(self.sphere_1_nc.accept(self.coords_100))
        self.assertFalse(self.sphere_20_nc.accept(self.coords_100))

    def test_positiveRadius(self):
        # test that using a negative value for radius correctly raises SignError
        self.assertRaises(exc.SignError, sphere, -1.0)
        self.assertRaises(exc.SignError, sphere, -1.0e-4)
        self.assertRaises(exc.SignError, sphere, -100.0)

    def test_numericalRadius(self):
        # test that using numbers for radius work, but that strings, etc. do not
        self.assertRaises(TypeError, sphere, "test_string")
        self.assertRaises(TypeError, sphere, "a")
        self.assertRaises(TypeError, sphere, "2test")
        self.assertRaises(TypeError, sphere, [2])
        self.assertRaises(TypeError, sphere, {"value": 2.0})
        self.assertIsInstance(sphere(2.0).radius2, type(1.0))
        self.assertIsInstance(sphere(2.0e-4).radius2, type(1.0))
        self.assertIsInstance(sphere(2).radius2, type(1.0))


if __name__ == '__main__':
    unittest.main()
