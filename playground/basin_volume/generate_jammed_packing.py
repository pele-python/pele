from __future__ import division
import numpy as np
import abc
import os
from pele.potentials import HS_WCA
from pele.optimize._quench import modifiedfire_cpp
from utils import *
import ConfigParser
import re

class _Generate_Jammed_Packing(object):
    """
    this is an abstract class that implements the basic components of a generate packing class,
    and declares a number of abstract methods which should be implemented in all inheriting classes
    *method to generate packing, this could be for example direct sampling, 
    sequential sampling,quench or LSA
    *nparticles: number of particles
    *bdim: dimensionality of the box
    *ndim: dimensionality of the problem (i.e. size of the coordinates array)
    *packing_frac: target jammed packing fraction
    *boxv: an array of size bdim that contains the vectors defining the box
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, packing_frac=0.65, packings_dir='packings'):
        self.packing_frac = packing_frac
        self.base_directory = os.path.join(os.getcwd(),'jammed_packings')
        self.packings_dir = os.path.join(os.getcwd(),packings_dir)
        self.configpath = os.path.join(packings_dir,'packings.config')
        self._import_packing_config_file()
        self.iteration = 0
        self.sca = -1
        
    def _import_packing_config_file(self):
        configf = ConfigParser.ConfigParser()
        configf.read(str(self.configpath))
        self.nparticles = configf.getint('PACKING','nparticles')
        self.bdim = configf.getint('PACKING','boxdim')
        assert(self.bdim==3) #currently PBC only implemented for 3d case
        self.ndim = self.nparticles * self.bdim
        boxv = configf.get('PACKING','boxv')
        self.boxv = np.array([float(x) for x in boxv.split()])
        self.imp_packing_frac = configf.getfloat('PACKING','packing_fraction')
    
    @abc.abstractmethod
    def initialise(self):
        """initialisation function"""
    
    @abc.abstractmethod
    def _import_packing_configuration(self, fname):
        """imports the coordinates and data relative to the shape of the particles
            this should be run in initialise()
        """
    @abc.abstractmethod
    def _generate_packing_coords(self):
        """function that generates the packing"""
    
    @abc.abstractmethod
    def _write_opengl_input(self, n):
        """writes a opengl input file, n is the unique identifier of the structure"""
        
    @abc.abstractmethod
    def _dump_configuration(self, n):
        """writes a configuration file, e.g .xyzd, n is the unique identifier of the structure"""
            
    def _print_initialise(self):
        base_directory = self.base_directory
        trymakedir(base_directory)
        self._print_parameters()
    
    def _print_parameters(self):
        """writes the simulation parameters"""
        fname = '{}/jammed_packings.config'.format(self.base_directory)
        f = open(fname,'w')
        f.write('#AUTOMATICALLY GENERATED FILE - DO NOT MODIFY BY HAND\n')
        f.write('#Generate_Jammed_Packings base class input parameters\n')
        f.write('[JAMMED_PACKING]\n')
        f.write('nparticles: {}\n'.format(self.nparticles))
        f.write('packing_fraction: {}\n'.format(self.packing_frac))
        f.write('boxdim: {}\n'.format(self.bdim))
        f.write('ndim: {}\n'.format(self.ndim))
        f.write('boxv: ')
        for val in self.boxv:
            f.write('{} '.format(val))
        f.write('\n')
        assert(self.sca >0)
        f.write('sca: {}\n'.format(self.sca))
        f.write('\n')
        f.close()
        
    def _print(self, n):
        """dump configuration and opengl input to packings directory
            n is the unique identifier of the structure
        """
        self._dump_configuration(n)
        self._write_opengl_input(n)
    
    @abc.abstractmethod
    def one_iteration(self,fname):
        """perform one iteration
        """
    def run(self):
        """run generate packings"""
        for fname in os.listdir(self.packings_dir):
            if 'xyzd' in fname:
                print fname
                self.one_iteration(fname)
            
class HS_Generate_Jammed_Packing(_Generate_Jammed_Packing):
    """
    *this class generates packings and identifies rattlers by computing the hessian eigenvalues for each particle
    *in the equilibrium jammed structure. A .xyzdr file is produced that contains the 3 system coordinates, the particle 
    * diameter and if not it's a rattler (0 if a rattler, 1 otherwise)
    *PARAMETERS
    *hs_radii: array with the radii of the particles, if none sample particle sizes from a normal distribution
    *mu: average particle size, passable to normal distribution
    *sig: standard deviaton of normal distribution from which to sample particles
    *sca: determines % by which the hs is inflated
    *eps: LJ interaction energy of WCA part of the HS potential
    """    
    def __init__(self, packing_frac=0.7, rattler_eval_tol=1.,packings_dir='packings'):
        super(HS_Generate_Jammed_Packing,self).__init__(packing_frac=packing_frac, packings_dir=packings_dir)
        
        ##constants#
        self.eps = 1.
        self.rattler_eval_tol = rattler_eval_tol 
        ############
    
    def initialise(self):
        self._compute_sca()
        self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxvec=self.boxv)
        self.rattlers = np.empty(self.nparticles,dtype='d')
        self._print_initialise()
    
    def one_iteration(self,fname):
        """perform one iteration
        """
        self._import_packing_configuration(fname)
        #initialise needs to import at least one configuration to compute sca
        if self.iteration is 0:
            self.initialise()
        self._generate_packing_coords()
        self._find_rattlers()
        #strips the integer unique identifier out of fname
        n = int(re.search(r'\d+',fname).group())
        self._print(n)
        self.iteration+=1
    
    def _find_rattlers(self):
        hess = self.potential.getHessian(self.coords)
        for i in xrange(self.nparticles):
            i1 = self.bdim*i
            hess_block = hess[i1:i1+self.bdim,i1:i1+self.bdim]
            w, v = np.linalg.eig(hess_block)
            rattler = np.less_equal(np.absolute(w),self.rattler_eval_tol)
            if True in rattler:
                self.rattlers[i] = 0.
                print 'zero eigenvalue, particle {}'.format(i)
                print w
            else:
                self.rattlers[i] = 1.
    
    def _generate_packing_coords(self):
        """quenches the imported structure using FIRE"""
        res = modifiedfire_cpp(self.coords,self.potential,nsteps=1e6, tol=1e-5)
        assert(res.success == True)
        self.coords = res.coords
        self.energy = res.energy
        no_overlap = self._check_overlaps()
        assert(no_overlap) #asserts that none of the hard sphere is overlapping
    
    def _get_particles_volume(self):
        """returns volume of n=self.bdim dimensional sphere"""
        volumes = 2*np.power(np.pi,self.bdim/2)*np.power(self.hs_radii,self.bdim)/(self.bdim*gamma(self.bdim/2))
        vtot = np.sum(volumes)
        return vtot
    
    def _import_packing_configuration(self, fname):
        path = os.path.join(self.packings_dir,fname)
        self.coords, self.hs_radii = read_xyzd(path)
    
    def _compute_sca(self):
        ##test##
        vol_part = self._get_particles_volume()
        vol_box = np.prod(self.boxv)
        phi = vol_part/vol_box #instanteneous pack frac
        assert(phi - self.imp_packing_frac < 1e-4)
        ##endtest##
        ###r_soft = r_hs*(1+sca)
        self.sca = np.power(self.packing_frac/self.imp_packing_frac,1./self.bdim) - 1
        
    def _check_overlaps(self):
        """check that no two particles are overlapping (using nearest image convention)"""
        no_overlap = True
        for i in xrange(self.nparticles):
            if no_overlap == True:
                for j in xrange(self.nparticles):
                    dij = 0
                    for k in xrange(self.bdim):
                        #use distances to nearest image convention
                        dij += np.square((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) -
                                          cround((self.coords[i*self.bdim+k] - self.coords[j*self.bdim+k]) / self.boxv[k]) * self.boxv[k])
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
            
    def _correct_coords(self):
        """this function returns the nearest images in the central box, useful for dumping the configurations"""
        coords = self.coords.copy()
        put_in_box(coords,self.boxv)
        return coords
    
    def _dump_configuration(self,n):
        """write coordinates to file .xyzdr"""
        coords = self._correct_coords()
        directory = self.base_directory
        fname = "{0}/jammed_packing{1}.xyzdr".format(directory,n)
        f = open(fname,'w')
        for i in xrange(self.nparticles):
            f.write('{:<12}\t{:<12}\t{:<12}\t{:<12}\t{:<12}\n'.format(coords[i*self.bdim],coords[i*self.bdim+1],
                                                               coords[i*self.bdim+2],self.hs_radii[i]*2,self.rattlers[i]))
        f.close()
    
    def _write_opengl_input(self,n):
        """write opengl input file"""
        coords = self._correct_coords()
        boxv = self.boxv
        colour = 14
        directory = self.base_directory
        fname = "{0}/jammed_packing{1}.dat".format(directory,n)
        f = open(fname,'w')
        f.write('{}\n'.format(self.nparticles))
        f.write('{} {} {}\n'.format(-boxv[0]/2,-boxv[1]/2,-boxv[2]/2))
        f.write('{} \t 0.0 \t 0.0\n'.format(boxv[0]))
        f.write('0.0 \t {} \t 0.0\n'.format(boxv[1]))
        f.write('0.0 \t 0.0 \t {}\n'.format(boxv[2]))
        for i in xrange(self.nparticles):
            for j in xrange(self.bdim):
                f.write('{}\t'.format(coords[i*self.bdim+j]))
            f.write('{}\t'.format(self.hs_radii[i]*2*(1.+self.sca)))
            f.write('{}\n'.format(colour-int(self.rattlers[i])))
        f.close()

            
if __name__ == "__main__":
    
    sim = HS_Generate_Jammed_Packing()
    sim.run()
    
    
        
                
            
              
                
                
                
