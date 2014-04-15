from __future__ import division
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.monte_carlo import Metropolis_MCrunner
from pele.parallel_tempering import MPI_PT_RLhandshake
from pele.potentials import Harmonic
        
if __name__ == "__main__":
    ndim = 3
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    Emax = 3
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    
    #Parallel Tempering
    mcrunner = Metropolis_MCrunner(potential, start_coords, niter=1e5, stepsize=1, hEmax = 100, adjustf = 0.9, 
                                   adjustf_niter = 5000, radius=100000)
    ptrunner = MPI_PT_RLhandshake(mcrunner, 0.2,1.6, max_ptiter=501, pfreq=10, verbose=False)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                
