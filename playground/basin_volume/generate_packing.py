from __future__ import division
import numpy as np
import abc
import os
from scipy.special import gamma
from playground.monte_carlo import Metropolis_MCrunner
from pele.potentials import HS_WCA, LJ
from pele.optimize._quench import lbfgs_cpp

#def volume_nball(radius,n):
#    volume = 2*np.power(np.pi,n/2)*np.power(radius,n)/(n*gamma(n/2))
#    return volume

def cround(r):
    if r > 0.0:
        r = np.floor(r + 0.5)
    else:
        r = np.ceil(r - 0.5)
    return r

def trymakedir(path):
    """this function deals with common race conditions"""
    while True:
        if not os.path.exists(path): 
            try:
                os.makedirs(path)
                break
            except OSError, e:
                if e.errno != 17:
                    raise
                # time.sleep might help here
                pass
        else:
            break

class _Generate_Packings(object):
    """
    *method to generate packing, this could be for example direct sampling, 
    sequential sampling,quench or LSA
    *nparticles: number of particles
    *bdim: dimensionality of the box
    *ndim: dimensionality of the problem (i.e. size of the coordinates array)
    *packing_frac: target packing fraction
    *boxv: an array of size bdim that contains the vectors defining the box
    *boxl: box side length
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, method, nparticles, bdim=3, boxv = None, boxl=1, packing_frac=0.4, max_iter = 1):
        self.method = method
        assert(bdim==3) #currently PBC only implemented for 3d case
        self.nparticles = nparticles
        self.bdim = bdim
        self.ndim = self.nparticles * self.bdim
        self.boxv = boxv
        self.boxl = boxl
        self.packing_frac = packing_frac
        self.initiliased=False
        self.base_directory = os.path.join(os.getcwd(),'packings')
        self.iteration = 0
        self.max_iter = max_iter
        self.box_resized = False
        self.initialised = False
        
    @abc.abstractmethod
    def initialise(self):
        """initialisation function"""
        
    @abc.abstractmethod
    def _get_particles_volume(self):
        """returns the total volume of the particles"""
        
    @abc.abstractmethod
    def _check_overlaps(self):
        """function to check that no two particles are overlapping"""
        
    @abc.abstractmethod
    def _generate_packing_coords(self):
        """function that generates the packing"""
    
    @abc.abstractmethod
    def _write_opengl_input(self):
        """writes a opengl input file"""
        
    @abc.abstractmethod
    def _dump_configuration(self):
        """writes a configuration file, e.g .xyzd"""
    
    def _resize_box(self):
        """adjust the box size to meet the target packing fraction"""
        vol_part = self._get_particles_volume()
        if (self.boxv == None):
            self.boxl = np.power(vol_part/self.packing_frac,1./3)
        else:
            vol_box = np.prod(self.boxv)
            phi = vol_part/vol_box #instanteneous pack frac
            a = np.power(phi/self.packing_frac,1./3)
            self.boxv *= a
        ###test###
        if (self.boxv == None):
            vol_box = self.boxl**self.bdim
        else:
            vol_box = np.prod(self.boxv)
        phi = vol_part/vol_box
        assert(phi - self.packing_frac < 1e-4)
        ##endtest## 
        self.box_resized = True
            
    def _print_initialise(self):
        base_directory = self.base_directory
        trymakedir(base_directory)
        self._print_parameters()
    
    def _print_parameters(self):
        """writes the simulation parameters"""
        fname = '{}/packings.config'.format(self.base_directory)
        f = open(fname,'w')
        f.write('#Generate_Packings base class input parameters\n')
        f.write('method: {}\n'.format(self.method))
        f.write('nparticles: {}\n'.format(self.nparticles))
        f.write('packing_fraction: {}\n'.format(self.packing_frac))
        f.write('boxdim: {}\n'.format(self.bdim))
        f.write('ndim: {}\n'.format(self.ndim))
        f.write('max_iter: {}\n'.format(self.max_iter))
        assert(self.box_resized)
        if self.boxv is None:
            f.write('boxl: {}\n'.format(self.boxl))
        else:
            f.write('boxv: ')
            for val in self.boxv:
                f.write('{} \t'.format(val))
            f.write('\n')
        f.close()
        
    def _print(self):
        """dump configuration and opengl input to packings directory"""
        self._dump_configuration()
        self._write_opengl_input()
    
    def one_iteration(self):
        """perform one iteration"""
        if self.initialised is not True:
            self.initialise()
        self._generate_packing_coords()
        self._print()
        self.iteration+=1
        
    def run(self):
        """run generate packings"""
        while self.iteration < self.max_iter:
            self.one_iteration()
            

class Generate_HS_Packing(_Generate_Packings):
    """
    *PARAMETERS
    *hs_radii: array with the radii of the particles, if none sample particle sizes from a normal distribution
    *mu: average particle size, passable to normal distribution
    *sig: standard deviaton of normal distribution from which to sample particles
    *sca: determines % by which the hs is inflated
    *eps: LJ interaction energy of WCA part of the HS potential, here irrelevant because 'sca' is set to 0
    """    
    def __init__(self, nparticles, method='quench', bdim=3, boxl=1, boxv=None, packing_frac=0.35, hs_radii=None, 
                 mu = 1, sig = 0.2, max_iter = 10):
        super(Generate_HS_Packing,self).__init__(method, nparticles, bdim=bdim, boxv = boxv, boxl=boxl,
                                                 packing_frac=packing_frac, max_iter = max_iter)
        ##constants#
        self.eps = 1.
        self.sca = 0.
        ############
        self.mu = mu
        self.sig = sig
        self.hs_radii = hs_radii
        self._sample_hs_radii()
        self._resize_box()
        
    def initialise(self):
        if self.method is 'quench':
            #this is necessary to initialise the radii if using the quench routine
            self._initialise_coords_quench()
        if self.boxv is None:
            self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxl=self.boxl)
        else:
            self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxv=self.boxv)
        self._print_initialise()
        self.initialised = True     
    
    def _get_particles_volume(self):
        """returns volume of n=self.bdim dimensional sphere"""
        volumes = 2*np.power(np.pi,self.bdim/2)*np.power(self.hs_radii,self.bdim)/(self.bdim*gamma(self.bdim/2))
        vtot = np.sum(volumes)
        return vtot
      
#    def _rescale_radii(self):
#        """rescale radii to meet target packing fraction"""
#        vol_box = np.power(self.boxl,self.bdim)
#        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
#        phi = vol_part/vol_box
#        self.hs_radii *= np.power(self.packing_frac/phi,1/self.bdim)
#        #test
#        vol_part = np.sum(volume_nball(self.hs_radii,self.bdim))
#        phi = vol_part/vol_box
#        assert(phi - self.packing_frac < 1e-4)
#        #endtest
    
    def _sample_hs_radii(self):
        if self.hs_radii is None:
            self.hs_radii = np.random.normal(self.mu,self.sig,self.nparticles)
        else:
            self.hs_radii = np.array(self.hs_radii,dtype='d')
        assert(self.hs_radii.all() > 0)
    
    def _check_overlaps(self):
        """check that no two particles are overlapping (using nearest image convention)"""
        if self.boxv is None:
            boxv = np.array([self.boxl for _ in xrange(self.bdim)])
        no_overlap = True
        for i in xrange(self.nparticles):
            if no_overlap == True:
                for j in xrange(self.nparticles):
                    dij = 0
                    for k in xrange(self.bdim):
                        #use distances to nearest image convention
                        dij += np.square((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) -
                                          cround((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) / boxv[k]) * boxv[k])
                    if i != j:
                        dij = np.sqrt(dij)
                        dmin = self.hs_radii[i]+self.hs_radii[j]
                        if dij - dmin <= 0:
                            print 'invalid configuration'
                            print 'atoms {} {} are overlapping'.format(i,j)
                            print 'real distance {}'.format(dij)
                            print 'min distance {}'.format(dmin)
                            no_overlap = False
                            break
            else:
                break
        return no_overlap
    
    def _sample_random_coords(self):
        """returns random coordinates for the particles uniformly distributed in the box"""
        if self.boxv is None:
            coords = np.random.rand(self.ndim)*self.boxl
        else:
            coords =  np.empty(self.ndim)
            for i in xrange(self.nparticles):
                for j in xrange(self.bdim):
                    coords[i*self.bdim+j] = np.random.rand()*self.boxv[j]
        return coords
    
    def _build_distance_matrix(self):
        distances = np.empty([self.nparticles,self.nparticles])
        for i in xrange(self.nparticles):
            for j in xrange(self.nparticles):
                dij = 0
                for k in xrange(self.bdim):
                    #use distances to closest image
                    dij += np.square((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) -
                                      cround((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) / self.boxl) * self.boxl)
                distances[i,j] = np.sqrt(dij)
        return distances
    
    def _generate_packing_coords(self):
        if self.method is 'quench':
            self._generate_packing_coords_quench()
        elif self.method is 'direct':
            self._generate_packing_coords_direct()
    
    def _generate_packing_coords_quench(self):
        """do a MCMC walk using the quenched coordinates"""
        if (self.iteration == 0):
            self.mcrunner = Metropolis_MCrunner(self.potential, self.coords, niter=1e5, stepsize=1e-4, adjustf = 0.9, 
                                       acceptance=0.2, adjustf_niter = 5000)
            self.energy = self.potential.getEnergy(self.coords)
        self.mcrunner.set_config(self.coords, self.energy)
        self.mcrunner.run()        
        self.coords, self.energy = self.mcrunner.get_config()
        
    def _initialise_coords_quench(self):
        """
        it generates an initial set of coordinates from a LJ quench,
        the LJ particles are then substitued by HS based on the size of
        the gap 
        """
        no_overlap = False
        if self.boxv is None:
            pot = LJ(sig=self.boxl,boxl=self.boxl)
        else:
            pot = LJ(sig=self.boxv[0],boxv=self.boxv) # choice of sigma might have to be different
        
        while no_overlap == False:
            no_overlap = True
            coords = self._sample_random_coords()
            res = lbfgs_cpp(coords,pot,nsteps=1000)
            #assert(res.success is True) #checks that a minimum configuration has been found
            self.coords = np.array(res.coords)
            print "generated new start coords "
            #build a matrix with the distances between particles i and j        
            distances = self._build_distance_matrix()
            #build an array with the weighted distance to neighbours, the shortest distance is 10 times heavier than the largest
            dmin = np.sort(distances,axis=1)
            
            if (self.nparticles > 12):
                neighbours = 12
            else:
                neighbours = self.nparticles
            
            CTE = np.exp( np.log(10) / (neighbours-1))
            weight = [CTE**i for i in xrange(neighbours)]
            weight = weight[::-1]
            weight.extend([0 for i in xrange(self.nparticles-neighbours)])
            #print 'weights',weight
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
            no_overlap = self._check_overlaps()
             
    def _generate_packing_coords_direct(self):
        """
        it generates an initial set of coordinates from a LJ quench,
        the LJ particles are then substitued by HS based on the size of
        the gap 
        """
        no_overlap = False
        while no_overlap == False:
            no_overlap = True 
            #coords = vector_random_uniform_hypersphere(self.ndim)#* np.sqrt(2*3) #coordinates sampled from Pow(ndim)
            self.coords = self._sample_random_coords()
            #assert(res.success is True) #checks that a minimum configuration has been found
            print "generated new start coords "
            #check that no two particles are overlapping (using nearest image convention)
            no_overlap = self._check_overlaps()
    
    def _dump_configuration(self):
        """write coordinates to file .xyzd"""
        directory = self.base_directory
        fname = "{0}/packing{1}.xyzd".format(directory,float(self.iteration))
        f = open(fname,'w')
        for i in xrange(self.nparticles):
            f.write('{:<12}\t{:<12}\t{:<12}\t{:<12}\n'.format(self.coords[i*self.bdim],self.coords[i*self.bdim+1],
                                                               self.coords[i*self.bdim+2],self.hs_radii[i]))
        f.close()
    
    def _write_opengl_input(self):
        """write open gl input file"""
        if self.boxv == None:
            boxv = [self.boxl for _ in xrange(self.bdim)]
        colour = 13
        directory = self.base_directory
        fname = "{0}/packing{1}.dat".format(directory,float(self.iteration))
        f = open(fname,'w')
        f.write('{}\n'.format(self.nparticles))
        f.write('0.0 0.0 0.0\n')
        f.write('{} \t 0.0 \t 0.0\n'.format(boxv[0]))
        f.write('0.0 \t {} \t 0.0\n'.format(boxv[1]))
        f.write('0.0 \t 0.0 \t {}\n'.format(boxv[2]))
        for i in xrange(self.nparticles):
            for j in xrange(self.bdim):
                f.write('{}\t'.format(self.coords[i*self.bdim+j]))
            f.write('{}\t'.format(self.hs_radii[i]*2))
            f.write('{}\n'.format(colour))
        f.close()

            
if __name__ == "__main__":
    
    nparticles = 30
    sim = Generate_HS_Packing(nparticles)
    sim.run()
#    sim.dump_coordinates_radii()
#    sim.write_opengl_input(0)
    #sim.run_mc()
    #sim.create_opengl_input(1)

    
    
        
                
            
              
                
                
                
