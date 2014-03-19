from __future__ import division
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from playground.monte_carlo import Metropolis_MCrunner
from playground.parallel_tempering import MPI_PT_RLhandshake
from pele.potentials import Harmonic
        
if __name__ == "__main__":
    ndim = 3
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    Emax = 3
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    
    #Parallel Tempering
    mcrunner = Metropolis_MCrunner(potential, start_coords, niter=1e6, stepsize=0.5, hEmax = 2000, adjustf = 0.95, 
                                   adjustf_niter = 5000, radius=1000)
    ptrunner = MPI_PT_RLhandshake(mcrunner, 0.2,1.6, max_ptiter=10, pfreq=2, verbose=True)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                