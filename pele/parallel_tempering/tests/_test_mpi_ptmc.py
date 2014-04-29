from __future__ import division
import numpy as np
import argparse
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.monte_carlo import Metropolis_MCrunner
from pele.parallel_tempering import MPI_PT_RLhandshake
from pele.potentials import Harmonic
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="do nested sampling on a Lennard Jones cluster")
    parser.add_argument("base_directory", type=str, help="directory in which to save results")
    #parser.add_argument("-K", "--nreplicas", type=int, help="number of replicas", default=300)
    args = parser.parse_args()
    
    ndim = 3
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    Emax = 3
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    path = args.base_directory
    
    #Parallel Tempering
    temperature=1.0
    stepsize=1
    niter=1e4
    mcrunner = Metropolis_MCrunner(potential, start_coords, temperature, stepsize, niter, hEmax = 100, adjustf = 0.9, 
                                   adjustf_niter = 3000, radius=100000)
    ptrunner = MPI_PT_RLhandshake(mcrunner, 0.2,1.6, max_ptiter=501, pfreq=10, base_directory=path, verbose=False)
    ptrunner.run()
            
            
            
            
            
    
    
    
        
                
            
              
                
                
                
