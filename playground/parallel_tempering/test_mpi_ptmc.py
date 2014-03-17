from __future__ import division
from mpi4py import MPI
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from playground.monte_carlo import Metropolis_MCrunner
from playground.parallel_tempering import MPI_PT_Simple
from pele.potentials import Harmonic, LJ
        
if __name__ == "__main__":
    #build harmonic potential
    #nparticles = 6
    #ndim = nparticles * 3
    ndim = 3
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    #potential = LJ()
    
    #build start configuration
    Emax = 3
    #start_coords = np.array([1,1,1,1,1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,1,-1])
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    
    #Parallel Tempering
    mcrunner = Metropolis_MCrunner(potential, start_coords, niter=1e5, stepsize=1, adjustf = 0.9, adjustf_niter = 1e5, radius=10000)
    ptrunner = MPI_PT_Simple(mcrunner, 0.2,1.6, max_ptiter=1000, pfreq=100)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                