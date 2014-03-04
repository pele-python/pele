from __future__ import division
import numpy as np
import random
from playground.monte_carlo import Metropolis_SHO_MCrunner
from playground.parallel_tempering import MPI_Parallel_Tempering
        
if __name__ == "__main__":
    mcrunner = Metropolis_SHO_MCrunner()
    ptrunner = MPI_Parallel_Tempering(mcrunner, 0.5,1.5)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                