import numpy as np
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic
from playground.monte_carlo import RandomCoordsDisplacement, MetropolisTest, RecordEnergyHistogram, MC

import abc
from pele.optimize import Result

class _base_mcrunner(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def get_results(self):
        """Must return a result object, generally must contain at least final configuration and energy"""

class test_mc(_base_mcrunner):
    
    def __init__(self, ndim = 20, niter=10000000, temperature=1.0,stepsize=0.8, binsize=0.01,k=1.0,Emin=7.0,Emax=11.0):
        self.niter = niter
        self.ndim = ndim
        self.temperature = temperature
        self.stepsize = stepsize/np.sqrt(ndim)
        self.binsize = binsize
        self.k = k
        self.Emax = Emax
        self.Emin = Emin
   
        #declare np array of desired dimensionality for origin of harmonic potential
        self.origin = np.zeros(self.ndim)
        #declare np array for actual coordinates
        self.coords = vector_random_uniform_hypersphere(self.ndim) * np.sqrt(2*self.Emax) #coordinates sampled from Pow(ndim)
        #declare np array where to save histogram at the end
        
        self.harmonic = Harmonic(self.origin,self.k)
        self.histogram = RecordEnergyHistogram(0,25,self.binsize)
        self.step = RandomCoordsDisplacement(self.ndim)
        self.metropolis = MetropolisTest()
        #set up mc
        self.mc = MC(self.harmonic, self.coords, temperature, self.stepsize)
        self.mc.set_takestep(self.step)
        self.mc.add_accept_test(self.metropolis)
        self.mc.add_action(self.histogram)
        self.result = Result()
        self.result.message = []
    
    def run(self):
        self.mc.run(self.niter)
    
    def get_results(self):
        """returns a results object"""
        res = self.result
        res.coords = self.mc.get_coords()
        res.energy = self.mc.get_energy()
        res.hist = self.histogram.get_histogram()
        res.accepted_frac = self.mc.get_accepted_fraction()
        return res
    
    def plot_histogram(self):
        import pylab as plt
        val = [i*self.binsize for i in xrange(len(self.result.hist))]
        plt.hist(val, weights=self.result.hist,bins=len(self.result.hist))
        plt.show()
    
    def print_dos_from_histogram(self):    
        """print the dos from histogram data"""
        
        # Open a file
        myfile = open("harmonic_dos.dat", "wb")
        
        E=self.Emin
        sumh = 0.0
        hist = self.result.hist
        n = int((self.Emax-self.Emin)/self.binsize)
        nstart = int(self.Emin/self.binsize)
        for i in xrange(nstart,nstart+n):
            sumh += hist[i] * np.exp(self.temperature*E) * self.binsize
            E += self.binsize
                      
        E = self.Emin
        for i in xrange(nstart,nstart+n):
            h = hist[i] * np.exp(self.temperature*E) / sumh
            myfile.write("{0} \t {1} \n".format(E,h))
            E += self.binsize
        
        myfile.close()
        
    def print_analytical_dos(self):
        """prints the analytical dos"""
        
        myfile = open("analytical_dos.dat", "wb")
        
        E=self.Emin
        n = int((self.Emax-self.Emin)/self.binsize)
        nstart = int(self.Emin/self.binsize)
        sumh=0.0
        for i in xrange(nstart,nstart+n):
            sumh += pow(E,self.ndim/2.0-1.0) * self.binsize
            E += self.binsize
              
        E = self.Emin
        for i in xrange(nstart,nstart+n):
            h = pow(E,self.ndim/2.0-1.0)/sumh
            myfile.write("{0} \t {1} \n".format(E,h))
            E += self.binsize
        
        myfile.close()

if __name__ == "__main__":
    test = test_mc()
    test.run()
    results = test.get_results()
    test.print_dos_from_histogram()
    test.print_analytical_dos()
    test.plot_histogram()

