from __future__ import absolute_import
import unittest
import numpy as np
import os

from pele.potentials import _morse_cpp
from pele.potentials.morse import Morse as PyMorse
from pele.utils.xyz import read_xyz
from . import _base_test

# class TestMorse(_base_test._BaseTest):
# def setUp(self):
# self.pot = PyMorse(rho=1.6047, r0=2.8970, A=0.7102)
# self.natoms = 13
# self.xrandom = np.random.uniform(-1,1,[3*self.natoms]) *5.

class TestMorse(_base_test._BaseTest):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_morse13_min.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin, rho, r0, A = list(map(float, xyz.title.split()[1::2]))

        self.natoms = self.xmin.size // 3
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        # self.pot = _morse_cpp.Morse(rho=rho, r0=r0, A=A)
        self.pot = PyMorse(rho=rho, r0=r0, A=A)


class TestMorse_CPP(_base_test._BaseTest):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        xyz = read_xyz(open(current_dir + "/_morse13_min.xyz", "r"))
        self.xmin = xyz.coords.reshape(-1).copy()
        self.Emin, rho, r0, A = list(map(float, xyz.title.split()[1::2]))

        self.natoms = self.xmin.size // 3
        self.xrandom = np.random.uniform(-1, 1, [3 * self.natoms]) * 5.
        self.pot = _morse_cpp.Morse(rho=rho, r0=r0, A=A)

        # self.pot = PyMorse(rho=rho, r0=r0, A=A)


def start_gui():
    from pele.systems import MorseCluster
    from pele.gui import run_gui

    natoms = 13
    system = MorseCluster(13, rho=1.6047, r0=2.8970, A=0.7102)
    run_gui(system)


if __name__ == "__main__":
    # start_gui()
    unittest.main()

