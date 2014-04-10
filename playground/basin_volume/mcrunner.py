import numpy as np
from pele.monte_carlo import _base_MCrunner, RandomCoordsDisplacement, MetropolisTest 
from pele.monte_carlo import CheckSameMinimum, AdjustStep, RecordDisp2Histogram
from pele.potentials import Harmonic, HS_WCA
from pele.optimize import ModifiedFireCPP

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

class HS_MCrunner(_base_MCrunner):
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
    def __init__(self, potential, coords, temperature=1.0, niter=1e5,
                  stepsize=1, hEmin=0, hEmax=100, hbinsize=0.01, radius=2.5,
                   acceptance=0., adjustf=0.9, adjustf_niter = 1e4, adjustf_navg = 100):
        #construct base class
        super(HS_MCrunner,self).__init__(potential, coords, temperature,
                                                  stepsize, niter)
                               
        #construct test/action classes       
        i32max = np.iinfo(np.int32).max
        
        self.adjust_step = AdjustStep(acceptance, adjustf, adjustf_niter, adjustf_navg)
        self.step = RandomCoordsDisplacement(self.ndim, np.random.randint(i32max))
        self.metropolis = MetropolisTest(np.random.randint(i32max))
        
        #set up pele:MC
        self.mc.set_takestep(self.step)
        self.mc.add_accept_test(self.metropolis)
        self.mc.add_action(self.adjust_step)
        
    def set_control(self, T):
        """set temperature, canonical control parameter"""
        self.temperature = T
        self.mc.set_temperature(T)
        
    
class BV_MCrunner(_base_MCrunner):
    """Basin Volume MCrunner
    *coords: initial coordinates, can be the same as origin
    *origin: jammed minimised structure
    *hs_radii: array of the radii of the particles
    *boxv: array with the box size lengths
    *rattlers: array of rattlers, if not rattler: 1 -> jammed dof
                                                  0 -> rattler dof
    *k: spring constant
    *temperature
    *niter: number of MC takesteps to perform
    *stepsize
    *Etol: tolerance with which a minimised structure is accepted
     when compared to origin energy
    *dtol: tolerance on the rms displacement of the minimised structure
     with respect to the origin coordinates    
    """
    def __init__(self, coords, origin, hs_radii, boxv, sca, rattlers=None, k=1.0, temperature=1.0, niter=1e5,
                  stepsize=0.01, dtol=1e-3, eps=1., hmin=0, hmax=100, hbinsize=0.01, acceptance=0.2, 
                  adjustf=0.9, adjustf_niter = 1e4, adjustf_navg = 100, opt_dtmax=1, opt_maxstep=0.5, opt_tol=1e-4,
                  opt_nsteps=1e5):
        #construct base class
        potential = Harmonic(origin,k,boxv)
        super(BV_MCrunner,self).__init__(potential, coords, temperature, stepsize, niter)
        
        self.origin = origin
        self.hs_radii = hs_radii
        self.boxv = boxv
        self.sca = sca
        self.k = k
        self.dtol = dtol
        self.eps = eps
        
        #manage array of rattlers, if not rattler: 1 -> jammed dof
        #                                          0 -> rattler dof 
        if (rattlers == None):
            self.rattlers = np.array([1. for _ in xrange(self.ndim)],dtype='d')
        else:
            self.rattlers = np.array(rattlers,dtype='d')
            assert(len(self.rattlers) == self.ndim)
            assert(self.rattlers.all() >= 0 and self.rattlers.all() <= 1)
        
        #construct gradient optimizer
        self.pot_optimizer = HS_WCA(self.eps, self.sca, self.hs_radii, boxvec=self.boxv)
        self.optimizer = ModifiedFireCPP(self.start_coords, self.pot_optimizer, dtmax=opt_dtmax, maxstep=opt_maxstep, 
                                         tol=opt_tol, nsteps=opt_nsteps)
                
        #construct test/action classes      
        i32max = np.iinfo(np.int32).max
        
        self.binsize = hbinsize
        self.histogram = RecordDisp2Histogram(self.origin, self.rattlers, hmin, hmax,self.binsize,adjustf_niter, self.boxv)
        self.conftest = CheckSameMinimum(self.optimizer, self.origin, self.hs_radii, self.boxv, self.rattlers, self.dtol)
        self.adjust_step = AdjustStep(acceptance, adjustf, adjustf_niter, adjustf_navg)
        self.step = RandomCoordsDisplacement(self.ndim, np.random.randint(i32max))
        self.metropolis = MetropolisTest(np.random.randint(i32max))
        
        #set up pele:MC
        self.mc.set_takestep(self.step)
        self.mc.add_accept_test(self.metropolis)
        self.mc.add_conf_test(self.conftest)
        self.mc.add_action(self.histogram)
        self.mc.add_action(self.adjust_step)
        
    def set_control(self, c):
        """set temperature, canonical control parameter"""
        self.k = c
        self.potential.set_k(c)
    
    def show_histogram(self):
        """shows the histogram"""
        import pylab as plt
        hist = self.histogram.get_histogram()
        val = [i*self.binsize for i in xrange(len(hist))]
        plt.hist(val, weights=hist,bins=len(hist))
        plt.show()
    
if __name__ == "__main__":
    #to run harmonic potential go to tests
    
    from pele.utils.rotations import vector_random_uniform_hypersphere
    from pele.optimize._quench import modifiedfire_cpp
    import time
    
    nparticles = 1
    ndim = nparticles * 3
    origin = np.array([0,0,0],dtype='d')
    #build start configuration
    Emax = 0.1
    start_coords = vector_random_uniform_hypersphere(ndim) * np.sqrt(2*Emax) #coordinates sampled from Pow(ndim)
    Harmonic(origin,1)
    res = modifiedfire_cpp(start_coords,Harmonic(origin,1))
    print res
    
#    print res.coords
    
    #Parallel Tempering
    test = BV_MCrunner(start_coords, origin, temperature=1, k=1, niter=1e5, hEmin=0,hEmax=100,
                       stepsize=0.5, adjustf = 0.9, adjustf_niter = 5000, radius=100)
    test.set_control(1)
    start=time.time()
    #test.run()
    end=time.time()
    print end-start
    test.show_histogram()
    