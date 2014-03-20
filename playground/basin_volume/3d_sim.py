from __future__ import division
import numpy as np
from scipy.special import gamma
from pele.utils.rotations import vector_random_uniform_hypersphere
from playground.monte_carlo import Metropolis_MCrunner
from playground.parallel_tempering import MPI_PT_RLhandshake
from pele.potentials import HS_WCA, WCA, LJ
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp

def volume_nball(radius,n):
    volume = 2*np.power(np.pi,n/2)*np.power(radius,n)/(n*gamma(n/2))
    return volume

class HSWCA_fluid_sim(object):
    """
    *PARAMETERS
    *nparticles
    *bdim is the dimensionality of the box
    *ndim is the dimensionality of the problem (i.e. size of the coordinates array)
    *packing fraction of hs
    *boxl: box side length
    *hs_radii: array with the radii of the particles, if none sample particle sizes from a normal distribution
    *sca: scaling factor that determines the thickness of the shell of the WCA+HS particles
    """    
    def __init__(self, nparticles, bdim=3, boxl=1, hs_radii=None, packing_frac=0.45, sca=0.1, eps=1.):
    
        self.nparticles = nparticles 
        self.bdim = bdim
        self.ndim = self.nparticles * self.bdim
        self.boxl = boxl
        self.packing_frac = packing_frac
        self.eps = eps
        self.sca = sca
        if hs_radii is None:
            self.hs_radii = np.absolute(np.random.randn(nparticles))
        else:
            self.hs_radii = hs_radii
    
    def initialise(self):
        self._rescale_radii()
        self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxl=self.boxl)
        self._generate_coords(self.potential)        
    
    def _rescale_radii(self):
        """rescale radii to meet target packing fraction"""
        vol_box = np.power(self.boxl,self.bdim)
        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
        phi = vol_part/vol_box
        self.hs_radii *= np.power(self.packing_frac/phi,1/self.bdim)
        #test
        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
        phi = vol_part/vol_box
        assert(phi - self.packing_frac < 1e-4)
        #endtest
    
    def _generate_coords(self,pot):
        """
        it generates an initial set of coordinate from which
        to start the simulation
        """
        coords = np.array([1.,1.,1.,1.,1.,0.])
        res = lbfgs_cpp(coords,pot)
        assert(res.success is True)
        self.coords = np.array(res.coords)
        dl = self.coords[:3] - self.coords[3:]
        d = np.sum(dl*dl)
        assert(np.sum(self.hs_radii)<d)
    
    def run_mc(self):
        mcrunner = Metropolis_MCrunner(self.potential, self.coords, niter=1e6, stepsize=0.01, adjustf = 0.9, 
                                       adjustf_niter = 5000)
        mcrunner.run()
        
        
#        Emax = 3
#        start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
#        
#        #Parallel Tempering
            
            
            
if __name__ == "__main__":
    
    nparticles = 2
    sim = HSWCA_fluid_sim(nparticles)
    sim.initialise()
    sim.run_mc()
        
                
            
              
                
                
                
