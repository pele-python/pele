from __future__ import division
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from playground.monte_carlo import Metropolis_MCrunner
from playground.parallel_tempering import MPI_PT_Simple
from pele.potentials import Harmonic
        
if __name__ == "__main__":
    #build harmonic potential
    ndim = 20
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    
    #build start configuration
    Emax = 25
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    
    #Parallel Tempering
    mcrunner = Metropolis_MCrunner(potential, start_coords)
    ptrunner = MPI_PT_Simple(mcrunner, 0.5,1.5)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                