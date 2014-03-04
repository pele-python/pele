from __future__ import division
import numpy as np
import random
from playground.monte_carlo import Metropolis_SHO_MCrunner
from playground.parallel_tempering import MPI_PT_Simple
        
if __name__ == "__main__":
    mcrunner = Metropolis_SHO_MCrunner()
    ptrunner = MPI_PT_Simple(mcrunner, 0.5,1.5)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                