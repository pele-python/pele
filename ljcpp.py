from math import *
import numpy as np
import ljcpp


class LJpshift:
    """binary lennard jones potential with smooth cutoff"""
    def __init__(self):
        # self.natoms = natoms

    def getEnergy(self, coords):
        print "getting energy only"
        E = ljcpp.energy(coords)
        return E

    def getEnergyGradient(self, coords):
        grad=np.zeros(coords.shape[0], np.float64)
        gE = ljcpp.energy(coords, grad)
        return E, V
