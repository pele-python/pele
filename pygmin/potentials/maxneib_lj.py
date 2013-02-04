"""
lj potential with the number of near neighbors restricted.
"""
import numpy as np

from pygmin.potentials import BasePotential
from pygmin.systems import LJCluster
import fortran.maxneib_lj as fortranpot



class MaxNeibsLJ(BasePotential):
    """
    atoms interact with lj potential, with an additional energy penalty for too many neighbors
    
    use the fermi function to get a continuous measure of whether two atoms are neighbors.  this gives
    a continuous measure of the number of neighbors of an atoms.
    
    use the fermi function again to calculate the energy penalty if an atom has too many neighbors.
    
    Parameters
    ----------
    eps : float
        lj eps
    sig : float
        lj sig
    boxl : float
        box length if periodic
    max_neibs : float 
        maximum number of allowed neighbors
    neib_crossover : float
        width of neighbor energy crossover function
    rneib : float
        distance to consider two atoms neighbors
    rneib_crossover : float
        width of rneib crossover region
    epsneibs : float
        energy scale of the neighbor penalty function
    
    
    """
    def __init__(self, eps=1.0, sig=1.0, boxl=None,
                 max_neibs=5.,
                 neib_crossover=.4,
                 rneib=1.4,
                 rneib_crossover=0.08,
                 epsneibs=5.,
                 ):
        """ simple lennard jones potential"""
        self.sig = sig
        self.eps = eps
        
        if boxl is None:
            self.boxl = 1e100
            self.periodic = False
        else:
            self.boxl = boxl
            self.periodic = True
        
        self.rneib = rneib
        self.rneib_crossover = rneib_crossover
        self.max_neibs = max_neibs
        self.neib_crossover = neib_crossover
        self.epsneibs = epsneibs

    def getEnergy(self, coords):
        E = fortranpot.maxneib_ljenergy(
                coords, self.eps, self.sig, self.periodic, self.boxl, 
                self.rneib, self.rneib_crossover, self.max_neibs, self.neib_crossover, 
                self.epsneibs)
        return E
    def getEnergyGradient(self, coords):
        E, grad = fortranpot.maxneib_ljenergy_gradient(
                coords, self.eps, self.sig, self.periodic, self.boxl, 
                self.rneib, self.rneib_crossover, self.max_neibs, self.neib_crossover, 
                self.epsneibs)
        return E, grad


class MaxNeibsLJSystem(LJCluster):
    def __init__(self, natoms, **potkwargs):
        super(MaxNeibsLJSystem, self).__init__(natoms)
        self.potkwargs = potkwargs
        self.params.gui.basinhopping_nsteps = 1000
    def __call__(self):
        return self
    
    def get_potential(self):
        return MaxNeibsLJ(**self.potkwargs)


def run_gui(system):
    import pygmin.gui.run as gr
    gr.run_gui(system)


if __name__ == "__main__":
    natoms = 20
    system = MaxNeibsLJSystem(natoms, max_neibs=3, rneib=1.7)
    
    coords = system.get_random_configuration()
    pot = system.get_potential()
    E = pot.getEnergy(coords)
    print "energy", E
    pot.test_potential(coords)
    if True:
        coords = system.get_random_minimized_configuration()[0]
        pot.test_potential(coords)
    
    run_gui(system)
