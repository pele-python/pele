import numpy as np
from playground.monte_carlo import _base_MCrunner, RandomCoordsDisplacement, MetropolisTest 
from playground.monte_carlo import CheckSphericalContainer, AdjustStep, RecordEnergyHistogram
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic

"""
pele::MCrunner

Specific implementations of MCrunners, generally they should follow this pattern:
* construct _base_MCrunner
* construct takestep, accept test, configuration test, action classes
* add these to the pele::MC class
* write a set_control function, for example you may want to set the temperature 
  (this is done this way to be compatible with the MPI replica exchange/parallel tempering
  implementation)
* add other functionalities that you may find desirable, e.g. dump histogram to file
"""

class Metropolis_MCrunner(_base_MCrunner):
    """This class is derived from the _base_MCrunner abstract
     method and performs Metropolis Monte Carlo. This particular implementation of the algorithm: 
     * runs niter steps per run call 
     * takes steps by sampling a random vector in a n dimensional hypersphere (n is the number of coordinates);
     * adjust the step size for the first adjustf_niter steps (averaging the acceptance for adjust_navg steps
       and adjusting the stepsize by a factor of 'adjustf') to meet some target acceptance 'acceptance'.
     * configuration test: accept if within a spherical box of radius 'radius'
     * acceptance test: metropolis for some particular temperature
     * record energy histogram (the energy histogram is resizable, but the bounds are defined by hEmin and hEmax,
       furthermore the bin size is set with hbinsize. Care must be taken because the array is resizable, if the step size
       is small and extremely high or low energies are sampled the memory for the histogram will be reallocated and this 
       might cause a badalloc error, if trying to allocate a huge array. If you are sampling unwanted extremely high or low energies
       then you might want to add a pele::EnergyWindow test that guarantees to keep you within a specific energy range and/or 
       make the stepsize larger or you might want to re-think about your simulation. Generally you shouldn't be 
       spanning energies that differ by several orders of magnitude, if that is the case, resizable or not resizable arrays are
       not the problem, you'd be incurring in memory issues no matter what you do, unless you write to disk at every iteration)
    """
    def __init__(self, potential, coords, temperature=1.0, niter=1e5,
                  stepsize=1, hEmin=0, hEmax=100, hbinsize=0.01, radius=2.5,
                   acceptance=0.5, adjustf=0.9, adjustf_niter = 1e4, adjustf_navg = 100):
        #construct base class
        super(Metropolis_MCrunner,self).__init__(potential, coords, temperature,
                                                  stepsize, niter)
                               
        #construct test/action classes       
        self.binsize = hbinsize
        self.histogram = RecordEnergyHistogram(hEmin,hEmax,self.binsize)
        self.adjust_step = AdjustStep(acceptance, adjustf, adjustf_niter, adjustf_navg)
        self.step = RandomCoordsDisplacement(self.ndim)
        self.metropolis = MetropolisTest()
        self.conftest = CheckSphericalContainer(radius)
        
        #set up pele:MC
        self.mc.set_takestep(self.step)
        self.mc.add_accept_test(self.metropolis)
        self.mc.add_conf_test(self.conftest)
        self.mc.add_action(self.histogram)
        self.mc.add_action(self.adjust_step)
        
    def set_control(self, T):
        """set temperature, canonical control parameter"""
        self.temperature = T
        self.mc.set_temperature(T)
    
    def dump_histogram(self, fname):
        """write histogram to fname"""
        Emin, Emax = self.histogram.get_Ebounds()
        histl = self.histogram.get_histogram()
        hist = np.array(histl)
        Energies, step = np.linspace(Emin,Emax,num=len(hist),endpoint=False,retstep=True)
        assert(abs(step - self.binsize) < self.binsize/100)
        np.savetxt(fname, np.column_stack((Energies,hist)), delimiter='\t')
    
    def get_histogram(self):
        """returns a energy list and a histogram list"""
        Emin, Emax = self.histogram.get_Ebounds()
        histl = self.histogram.get_histogram()
        hist = np.array(histl)
        Energies, step = np.linspace(Emin,Emax,num=len(hist),endpoint=False,retstep=True)
        assert(abs(step - self.binsize) < self.binsize/100)
        return Energies, hist
        
    def show_histogram(self):
        """shows the histogram"""
        import pylab as plt
        hist = self.histogram.get_histogram()
        val = [i*self.binsize for i in xrange(len(hist))]
        plt.hist(val, weights=hist,bins=len(hist))
        plt.show()
    
#    def print_dos_from_histogram(self,Emin=7.0,Emax=11.0):    
#        """print the dos from histogram data"""
#        
#        # Open a file
#        myfile = open("harmonic_dos.dat", "wb")
#        
#        E=Emin
#        sumh = 0.0
#        histl = self.histogram.get_histogram()
#        hist = np.array(histl)
#        n = int((Emax-Emin)/self.binsize)
#        nstart = int(Emin/self.binsize)
#        for i in xrange(nstart,nstart+n):
#            sumh += hist[i] * np.exp(self.temperature*E) * self.binsize
#            E += self.binsize
#                      
#        E = Emin
#        for i in xrange(nstart,nstart+n):
#            h = hist[i] * np.exp(self.temperature*E) / sumh
#            myfile.write("{0} \t {1} \n".format(E,h))
#            E += self.binsize
#        
#        myfile.close()
#        
#    def print_analytical_dos(self,Emin=7.0,Emax=11.0):
#        """prints the analytical dos"""
#        
#        myfile = open("analytical_dos.dat", "wb")
#        
#        E=Emin
#        n = int((Emax-Emin)/self.binsize)
#        nstart = int(Emin/self.binsize)
#        sumh=0.0
#        for i in xrange(nstart,nstart+n):
#            sumh += pow(E,self.ndim/2.0-1.0) * self.binsize
#            E += self.binsize
#              
#        E = Emin
#        for i in xrange(nstart,nstart+n):
#            h = pow(E,self.ndim/2.0-1.0)/sumh
#            myfile.write("{0} \t {1} \n".format(E,h))
#            E += self.binsize
#        
#        myfile.close()
        
if __name__ == "__main__":
    #build harmonic potential
    ndim = 3
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    
    cvs = []
    #build start configuration
    temperatures = [0.2,0.27,0.362,0.487,0.65,0.88,1.18,1.6]
    stepsizes=[0.378,0.467,0.467,0.577,0.712,0.88,0.88,1.34]
    for k,T in enumerate(temperatures):
        for j in xrange(100):
            Emax = 3
            start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
            
            #MCMC 
            test = Metropolis_MCrunner(potential, start_coords,  temperature=T, niter=1e7, stepsize=stepsizes[k], adjustf = 0.9, adjustf_niter = 1000, radius=10000)
            test.run()
            #collect the results
            #test.show_histogram()
            binenergy, hist = test.get_histogram()
            
            average = np.average(binenergy,weights=hist)
                
            average2 = np.average(np.square(binenergy),weights=hist)
                
            cv =  (average2 - average**2)/(T**2) + 1.5
            
            print 'heat capacity ',cv
    
            cvs.append(cv)
        temp = np.empty(len(cvs))
        temp.fill(T)
        np.savetxt('heat_capacities_{}'.format(T), np.column_stack((temp,cvs)),delimiter='\t')
        cvs =[]