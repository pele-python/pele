import numpy as np
from pele.monte_carlo import _BaseMCRunner, RandomCoordsDisplacement, MetropolisTest 
from pele.monte_carlo import CheckSphericalContainer, AdjustStep, RecordEnergyHistogram

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

class Metropolis_MCrunner(_BaseMCRunner):
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
     * NOTE: some of the modules (e.g. take step and acceptance tests) require to be seeded. Users are free to do this as they think
     * is best, here we generate a random integer in [0,i32max) where i32max is the largest signed integer, for each seed. Each module
     * has a separate rng engine, therefore it's best if each receives a different randomly sampled seed
    """
    def __init__(self, potential, coords, temperature, stepsize, niter, 
                 hEmin=0, hEmax=100, hbinsize=0.01, radius=2.5,
                 acceptance=0.5, adjustf=0.9, adjustf_niter = 1e4, adjustf_navg = 100):
        #construct base class
        super(Metropolis_MCrunner,self).__init__(potential, coords, temperature,
                                                  stepsize, niter)
                               
        #construct test/action classes       
        i32max = np.iinfo(np.int32).max
        
        self.binsize = hbinsize
        self.histogram = RecordEnergyHistogram(hEmin,hEmax,self.binsize, adjustf_niter)
        self.adjust_step = AdjustStep(acceptance, adjustf, adjustf_niter, adjustf_navg)
        self.step = RandomCoordsDisplacement(self.ndim, np.random.randint(i32max))
        self.metropolis = MetropolisTest(np.random.randint(i32max))
        self.conftest = CheckSphericalContainer(radius)
        
        #set up pele:MC
        self.set_takestep(self.step)
        self.add_accept_test(self.metropolis)
        self.add_conf_test(self.conftest)
        self.add_action(self.histogram)
        self.add_action(self.adjust_step)
        
    def set_control(self, T):
        """set temperature, canonical control parameter"""
        self.temperature = T
        self.set_temperature(T)
    
    def dump_histogram(self, fname):
        """write histogram to fname"""
        Emin, Emax = self.histogram.get_bounds_val()
        histl = self.histogram.get_histogram()
        hist = np.array(histl)
        Energies, step = np.linspace(Emin,Emax,num=len(hist),endpoint=False,retstep=True)
        assert(abs(step - self.binsize) < self.binsize/100)
        np.savetxt(fname, np.column_stack((Energies,hist)), delimiter='\t')
    
    def get_histogram(self):
        """returns a energy list and a histogram list"""
        Emin, Emax = self.histogram.get_bounds_val()
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
    #to run harmonic potential go to tests
    
    from pele.optimize import LBFGS_CPP
    from pele.potentials import LJ
    import time
    
    nparticles = 31
    ndim = nparticles * 3
    pot = LJ()
    #build start configuration
#    Emax = 3
#    xrand = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
#    e, grad = pot.getEnergyGradient(xrand)
#    lbfgs = LBFGS_CPP(xrand, pot, energy=e, gradient=grad)
#    res = lbfgs.run()
#    print res.coords
    #start_coords6 = np.array([1,1,1,1,1,-1,1,-1,1,-1,1,1,1,-1,-1,-1,1,-1])
    #lj31 stable coords, from minimisation, convenient keeping it in this form for profiling
    start_coords = np.array([-0.77948249,-0.1563986,1.14725592,1.13438658,0.68768745,-1.01400641,
                               0.03889983,0.2795145,0.50005794,1.55312878,-0.80384468,0.0632129,
                               1.28113216,1.08082843,0.01281633,-0.51859703,-1.14310736,-0.642079,
                               -0.74400707,0.77214221,-0.78159587,-0.39137658,0.88003139,1.31292537,
                               -1.41406151,0.73317798,0.91291251,-1.36912846,-0.86828268,0.03339192,
                               -0.28968205,-0.50094652,-1.53635188,0.68474233,-0.60974182,0.72916478,
                               1.06365511,0.41981176,0.88587376,0.52582521,-0.97592561,-0.30476306,
                               -0.11140428,-0.10604705,-0.53162998,-0.32339665,-0.71798417,0.33592738,
                               0.3464601,-1.37722735,-1.30968712,-0.12150755,-1.02059093,1.37780291,
                               -1.15472311,-0.25677839,-0.88824791,0.25781619,0.90928202,-0.35987022,
                               0.45434582,1.56354974,-1.20590082,0.41425323,1.2978406,0.67123977,
                               -0.96365292,0.14063053,0.1074299,0.11754687,0.52112844,-1.41019987,
                               0.89892003,0.05030545,-0.14158607,-0.61370579,1.16407401,0.25903675,
                               0.76311679,-0.34748748,-1.16950066,0.2565447,0.01168485,1.54106553,
                               0.24769915,-1.6362588,0.5225262,-0.85826368,-1.5714381,0.73425918,
                               -1.64261439,0.95488524,-0.16297967]) 
    
    #Parallel Tempering
    temperature=0.2
    niter=1e7
    stepsize=0.5
    test = Metropolis_MCrunner(pot, start_coords,  temperature, stepsize, niter, 
                               hEmin=-140, adjustf = 0.9, adjustf_niter = 5000, radius=3)
    start=time.time()
    test.run()
    end=time.time()
    print end-start
    
    