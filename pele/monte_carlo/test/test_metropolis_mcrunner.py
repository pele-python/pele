from __future__ import division
import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic
from pele.monte_carlo import Metropolis_MCrunner
import unittest
import logging


class TestMetropolis(unittest.TestCase):
    
    def test_heat_capacity(self):
        self.ndim = 3
        self.k=1
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin,self.k)
        self.Emax = 2 #this choice is fundamentally arbitrary, it's only used to generate the initial configuration
        self.start_coords = vector_random_uniform_hypersphere(self.ndim) * np.sqrt(2*self.Emax) #coordinates sampled from Pow(ndim)
               
        temperatures = [0.2,0.27,0.362,0.487,0.65,0.88,1.18,1.6]
                
        for T in temperatures:
            temperature=T
            stepsize=0.5
            niter=1e6
            mcrunner = Metropolis_MCrunner(self.potential, self.start_coords, temperature, stepsize, niter, hEmax = 100, adjustf = 0.9, 
                                           adjustf_niter = 1000, radius=100000)
            #MCMC 
            mcrunner.run()
            
            #collect the results
            binenergy, hist = mcrunner.get_histogram()
            
            average = np.average(binenergy,weights=hist)
                
            average2 = np.average(np.square(binenergy),weights=hist)
                
            cv =  (average2 - average**2)/(T**2) + self.ndim/2
            
            self.assertLess(cv-self.ndim,0.1,'failed for temperature {}, cv = {}'.format(T,cv))

if __name__ == "__main__":
    logging.basicConfig(filename='Metropolis_mcrunner.log',level=logging.DEBUG)
    unittest.main()
