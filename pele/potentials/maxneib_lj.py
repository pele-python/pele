"""
lj potential with the number of near neighbors restricted.
"""
import numpy as np

from pele.potentials import BasePotential
from pele.systems import LJCluster
import fortran.maxneib_lj as fortranpot

__all__ = ["MaxNeibsLJ", "MaxNeibsLJSystem"]


class MaxNeibsLJ(BasePotential):
    """
    atoms interact with lj potential, with an additional energy penalty for too many neighbors
    
    Notes
    -----
    in addition to the LJ energy the penalty term is::
        
        E = sum_i F(ni, max_neibs, neib_crossover)
        
    where `F` is the fermi function defined as::
    
        F(x, mu, T) = 1. / (exp(-(x-mu)/T) + 1.)
    
    and `ni` is a continuous measure of the number of neighbors of the i'th
    atom::
    
        ni = sum_j F(rij, rneib, rneib_crossover)
    
    where `rij` is the separation of atoms `i` and `j`.  The other parameters
    above are defined in the Parameters section

    The Fermi function is used in order to get a continuous measure of whether
    two atoms are neighbors.  It is also used again to calculate the energy
    penalty if an atom has too many neighbors.
    
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
    
    
    See Also
    --------
    MaxNeibsBLJ
    
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
        if self.periodic:
            coords -= np.round(coords / self.boxl) * self.boxl
        E, grad = fortranpot.maxneib_ljenergy_gradient(
                coords, self.eps, self.sig, self.periodic, self.boxl, 
                self.rneib, self.rneib_crossover, self.max_neibs, self.neib_crossover, 
                self.epsneibs)
        return E, grad


class MaxNeibsLJSystem(LJCluster):
    def __init__(self, natoms, **potkwargs):
        super(MaxNeibsLJSystem, self).__init__(natoms)
        self.potkwargs = potkwargs
        self.params.gui.basinhopping_nsteps = 300
    def __call__(self):
        return self
    
    def get_potential(self):
        return MaxNeibsLJ(**self.potkwargs)


def run_gui(system):
    import pele.gui.run as gr
    gr.run_gui(system)


if __name__ == "__main__":
    periodic = True
    natoms = 20
    if periodic:
        rho = 0.5
        boxl = (float(natoms) / rho)**(1./3)
        print boxl
#        exit()
    else:
        boxl = None
    system = MaxNeibsLJSystem(natoms, max_neibs=3, rneib=1.7, boxl=boxl, epsneibs=12.)
    
    coords = system.get_random_configuration()
    pot = system.get_potential()
    E = pot.getEnergy(coords)
    print "energy", E
    pot.test_potential(coords)
    if True:
        coords = system.get_random_minimized_configuration()[0]
        pot.test_potential(coords)
#    exit(10)
    
    run_gui(system)
