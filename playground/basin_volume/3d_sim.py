from __future__ import division
import numpy as np
import heapq
from scipy.special import gamma
from pele.utils.rotations import vector_random_uniform_hypersphere
from playground.monte_carlo import Metropolis_MCrunner
from playground.parallel_tempering import MPI_PT_RLhandshake
from pele.potentials import HS_WCA, WCA, LJ
from pele.optimize._quench import lbfgs_cpp, modifiedfire_cpp

def volume_nball(radius,n):
    volume = 2*np.power(np.pi,n/2)*np.power(radius,n)/(n*gamma(n/2))
    return volume

def cround(r):
    if r > 0.0:
        r = np.floor(r + 0.5)
    else:
        r = np.ceil(r - 0.5)
    return r

class HSWCA_fluid_sim(object):
    """
    *PARAMETERS
    *nparticles
    *bdim is the dimensionality of the box
    *ndim is the dimensionality of the problem (i.e. size of the coordinates array)
    *packing fraction of hs
    *boxl: box side length
    *hs_radii: array with the radii of the particles, if none sample particle sizes from a normal distribution
    *sca: scaling factor that determines the thickness of the shell of the WCA+HS particles
    """    
    def __init__(self, nparticles, bdim=3, boxl=1, hs_radii=None, packing_frac=0.45, sca=0, eps=1.):
    
        assert(bdim==3) #currently PBC only implemented for 3d case
        self.nparticles = nparticles 
        self.bdim = bdim
        self.ndim = self.nparticles * self.bdim
        self.boxl = boxl
        self.packing_frac = packing_frac
        self.eps = eps
        self.sca = sca
        if hs_radii is None:
            self.hs_radii = np.random.randn(nparticles)+10 #distribution is shifted by 10
        else:
            self.hs_radii = hs_radii
        self.initialised = False
    
    def initialise(self):
        self._rescale_radii()
        self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxl=self.boxl)
        self._generate_coords(self.potential)
        self.initialised = True        
    
    def _rescale_radii(self):
        """rescale radii to meet target packing fraction"""
        vol_box = np.power(self.boxl,self.bdim)
        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
        phi = vol_part/vol_box
        self.hs_radii *= np.power(self.packing_frac/phi,1/self.bdim)
        #test
        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
        phi = vol_part/vol_box
        assert(phi - self.packing_frac < 1e-4)
        #endtest
    
    def _generate_coords(self,pot):
        """
        it generates an initial set of coordinate from which
        to start the simulation
        """
        no_overlap = False
        while no_overlap == False:
            no_overlap = True
            coords = vector_random_uniform_hypersphere(self.ndim) #* np.sqrt(2*3) #coordinates sampled from Pow(ndim)
            pot = LJ(sig=1,boxl=self.boxl)
            res = lbfgs_cpp(coords,pot,nsteps=1000)
            #assert(res.success is True) #checks that a minimum configuration has been found
            self.coords = np.array(res.coords)
            print "generated new start coords "
#            #build a matrix with the distances between particles i and j        
            distances = np.empty([self.nparticles,self.nparticles])
            for i in xrange(self.nparticles):
                for j in xrange(self.nparticles):
                    dij = 0
                    for k in xrange(self.bdim):
                        #use distances to closest image
                        dij += np.square((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) -
                                          cround((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) / self.boxl) * self.boxl)
                    distances[i,j] = np.sqrt(dij)
            #build an array with the weighted distance to neighbours, the shortest distance is 10 times heavier than the largest
            #print distances
            dmin = np.sort(distances,axis=1)
            CTE = np.exp( np.log(10) / (self.nparticles-1))
            weight = [CTE**i for i in range(self.nparticles)]
            weight = weight[::-1]
            for i in xrange(6,len(weight)):
                weight[i] = 0
            print 'weights',weight
            dmin = np.average(dmin,axis=1,weights=weight)
            #sort and return a map of indices in descending order
            dmap = np.argsort(dmin)[::-1]
            #order particle sizes so that they are associated to coordinates with appropriate gaps
            #print 'old radii',self.hs_radii
            hs_radii = np.zeros(self.nparticles)
            sorted_radii = np.sort(self.hs_radii)[::-1]
            for i in xrange(self.nparticles):
                hs_radii[dmap[i]] = sorted_radii[i]
            self.hs_radii = hs_radii.copy()
            #print 'new radii',self.hs_radii
            
            #check that no two particles are overlapping (using nearest image convention)
            for i in xrange(self.nparticles):
                if no_overlap == True:
                    for j in xrange(self.nparticles):
                        dij = 0
                        for k in xrange(self.bdim):
                            #use distances to nearest image convention
                            dij += np.square((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) -
                                              cround((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) / self.boxl) * self.boxl)
                        if i != j:
                            dij = np.sqrt(dij)
                            dmin = self.hs_radii[i]+self.hs_radii[j]
                            if dij - dmin < 0:
                                print '{} {}'.format(i,j)
                                print 'real distance {}'.format(dij)
                                print 'min distance {}'.format(dmin)
                                no_overlap = False
                                break
                else:
                    break
    
    def run_mc(self):
        if self.initialised is not True:
            self.initialise()
        #res = modifiedfire_cpp(self.coords,self.potential,nsteps=1000)
        #self.coords = res.coords
        mcrunner = Metropolis_MCrunner(self.potential, self.coords, niter=1e5, stepsize=1e-4, adjustf = 0.9, 
                                       adjustf_niter = 5000)
        mcrunner.run()
        print 'stepsize',mcrunner.get_stepsize()
        print 'accepted fraction',mcrunner.get_accepted_fraction()
        print 'config',mcrunner.get_config()[0] - self.coords
        
        
            
            
            
if __name__ == "__main__":
    
    nparticles = 100
    sim = HSWCA_fluid_sim(nparticles)
    sim.initialise()
    sim.run_mc()

    
    
        
                
            
              
                
                
                
