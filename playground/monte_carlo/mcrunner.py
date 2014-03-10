import numpy as np
import abc
from pele.utils.rotations import vector_random_uniform_hypersphere
from pele.potentials import Harmonic
from playground.monte_carlo import RandomCoordsDisplacement, MetropolisTest, CheckSphericalContainer, AdjustStep, RecordEnergyHistogram, MC
from pele.optimize import Result

class _base_MCrunner(object):
    """
    Abstract method for MC runners, all MC runners should derive from this base class
    *potential should be constructed outside of this class and passed
    *coords are the initial coordinates
    *niter is the total number of MC iterations
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, potential, coords, temperature, stepsize, niter):
        self.potential = potential
        self.ndim = len(coords)
        self.start_coords = coords
        self.temperature = temperature
        self.stepsize = stepsize/np.sqrt(self.ndim)
        self.niter = niter
        self.mc = MC(self.potential, self.start_coords, self.temperature, self.stepsize)    
        self.result = Result()
        self.result.message = []
        
    def get_config(self):
        """Return the coordinates of the current configuration and its associated energy"""
        coords = self.mc.get_coords()
        energy = self.mc.get_energy()
        return coords, energy
    
    def set_config(self, coords, energy):
        """set current configuration and its energy"""
        self.mc.set_coordinates(coords, energy)
    
    def get_results(self):
        """Must return a result object, generally must contain at least final configuration and energy"""
        res = self.result
        res.coords = self.mc.get_coords()
        res.energy = self.mc.get_energy()
        res.accepted_frac = self.mc.get_accepted_fraction()
        return res
    
    def get_iterations_count(self):
        n = self.mc.get_iterations_count()
        return n
    
    def get_accepted_fraction(self):
        n = self.mc.get_accepted_fraction()
        return n
    
    @abc.abstractmethod
    def set_control(self):
        """set control parameter, this could be temperature or some other control parameter like stiffness of the harmonic potential"""
    
    def run(self):
        """run MCMC walk"""
        self.mc.run(self.niter)
        
class Metropolis_MCrunner(_base_MCrunner):
    """This class is derived from the _base_MCrunner abstract method and performs Metropolis Monte Carlo
    """
    def __init__(self, potential, coords, temperature=1.0, niter=1000000, stepsize=0.8, hEmin=0, hEmax=100, hbinsize=0.01, radius=5.0, acceptance=0.5, adjustf=0.9999, adjustf_niter = 5000):
        super(Metropolis_MCrunner,self).__init__(potential, coords, temperature, stepsize, niter)
                               
        #construct test/action classes       
        self.binsize = hbinsize
        self.histogram = RecordEnergyHistogram(hEmin,hEmax,self.binsize)
        self.adjust_step = AdjustStep(acceptance, adjustf, adjustf_niter)
        self.step = RandomCoordsDisplacement(self.ndim)
        self.metropolis = MetropolisTest()
        self.conftest = CheckSphericalContainer(radius)
        #set up mc
        self.mc.set_takestep(self.step)
        self.mc.add_accept_test(self.metropolis)
        self.mc.add_conf_test(self.conftest)
        self.mc.add_action(self.histogram)
        self.mc.add_action(self.adjust_step)
        
    def get_results(self):
        """returns a results object"""
        res = self.result
        res.coords = self.mc.get_coords()
        res.energy = self.mc.get_energy()
        res.accepted_frac = self.mc.get_accepted_fraction()
        res.hist = self.histogram.get_histogram()
        return res
    
    def set_control(self, T):
        self.temperature = T
        self.mc.set_temperature(T)
    
    def dump_histogram(self, fname):
        Emin, Emax = self.histogram.get_Ebounds()
        hist = np.array(self.result.hist)
        Energies, step = np.linspace(Emin,Emax,num=len(hist),endpoint=False,retstep=True)
        assert(abs(step - self.binsize) < self.binsize/100)
        #print "energies len {}, hist len {}".format(len(Energies),len(hist))
        np.savetxt(fname, np.column_stack((Energies,hist)), delimiter='\t') 
        
    def plot_histogram(self):
        import pylab as plt
        val = [i*self.binsize for i in xrange(len(self.result.hist))]
        plt.hist(val, weights=self.result.hist,bins=len(self.result.hist))
        plt.show()
    
    def print_dos_from_histogram(self,Emin=7.0,Emax=11.0):    
        """print the dos from histogram data"""
        
        # Open a file
        myfile = open("harmonic_dos.dat", "wb")
        
        E=Emin
        sumh = 0.0
        hist = self.result.hist
        n = int((Emax-Emin)/self.binsize)
        nstart = int(Emin/self.binsize)
        for i in xrange(nstart,nstart+n):
            sumh += hist[i] * np.exp(self.temperature*E) * self.binsize
            E += self.binsize
                      
        E = Emin
        for i in xrange(nstart,nstart+n):
            h = hist[i] * np.exp(self.temperature*E) / sumh
            myfile.write("{0} \t {1} \n".format(E,h))
            E += self.binsize
        
        myfile.close()
        
    def print_analytical_dos(self,Emin=7.0,Emax=11.0):
        """prints the analytical dos"""
        
        myfile = open("analytical_dos.dat", "wb")
        
        E=Emin
        n = int((Emax-Emin)/self.binsize)
        nstart = int(Emin/self.binsize)
        sumh=0.0
        for i in xrange(nstart,nstart+n):
            sumh += pow(E,self.ndim/2.0-1.0) * self.binsize
            E += self.binsize
              
        E = Emin
        for i in xrange(nstart,nstart+n):
            h = pow(E,self.ndim/2.0-1.0)/sumh
            myfile.write("{0} \t {1} \n".format(E,h))
            E += self.binsize
        
        myfile.close()
        
if __name__ == "__main__":
    #build harmonic potential
    ndim = 20
    k=1
    origin = np.zeros(ndim)
    potential = Harmonic(origin,k)
    
    #build start configuration
    Emax = 1
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    
    #MCMC 
    test = Metropolis_MCrunner(potential, start_coords)
    test.run()
    results = test.get_results()
    #test.print_dos_from_histogram()
    #test.print_analytical_dos()
    test.plot_histogram()