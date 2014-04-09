from __future__ import division
import numpy as np
import abc
import os
from pele.potentials import HS_WCA
from utils import *
import ConfigParser
import time
import re
import pylab

class _analyse_jammed_packings(object):
    """
    this is an abstract class that implements the basic components of a configure bv_mcrunner class,
    and declares a number of abstract methods which should be implemented in all inheriting classes
    *nparticles: number of particles
    *bdim: dimensionality of the box
    *ndim: dimensionality of the problem (i.e. size of the coordinates array)
    *packing_frac: target jammed packing fraction
    *boxv: an array of size bdim that contains the vectors defining the box
    """
        
    def __init__(self, packings_dir='jammed_packings'):
        self.base_directory = os.path.join(os.getcwd(),'analyse_jammed_packings')
        self.packings_dir = os.path.join(os.getcwd(),packings_dir)
        self.configpath = os.path.join(packings_dir,'jammed_packings.config')
        self._import_packing_config_file()
        self.eps=1.
        self.iteration = 0
        self.eigenvalues = []
        self.nbins = 100
        #construct mcrunner
        #self.coords is origin, set initial configuration and origin to be the same
        
    def _initialise(self):
        """initialisation function"""
        self._import_packing_config_file()
        #change directory only at the end of initialise
        self._print_initialise()
        os.chdir(self.base_directory)
        
    def _import_packing_config_file(self):
        configf = ConfigParser.ConfigParser()
        configf.read(str(self.configpath))
        self.nparticles = configf.getint('JAMMED_PACKING','nparticles')
        self.bdim = configf.getint('JAMMED_PACKING','boxdim')
        assert(self.bdim==3) #currently PBC only implemented for 3d case
        self.ndim = self.nparticles * self.bdim
        boxv = configf.get('JAMMED_PACKING','boxv')
        self.boxv = np.array([float(x) for x in boxv.split()])
        self.imp_packing_frac = configf.getfloat('JAMMED_PACKING','packing_fraction')
        self.sca = configf.getfloat('JAMMED_PACKING','sca')
        
    def _import_packing_configuration(self,fname):
        """imports the coordinates and data relative to the shinitape of the particles
            this should be run in initialise()
        """
        path = os.path.join(self.packings_dir,fname)
        self.coords, self.hs_radii = read_xyzd(path)
    
    def _histogram_eigenvalues(self):
        #self.eigenvalues = np.array(self.eigenvalues,dtype='d')
        self.histogram, bins = np.histogram(self.eigenvalues,bins=self.nbins)
        #print self.histogram
        #n, bins, patches = pylab.hist(self.eigenvalues,self.nbins,histtype='bar')
        width = bins[1] - bins[0]
        center = (bins[:-1] + bins[1:]) / 2
        pylab.bar(center, self.histogram, align='center', width=width)
        pylab.show()
    
    def one_iteration(self,fname):
        """compute hessian and its eigenvalues
        """
        self.hess_block = np.zeros((self.bdim,self.bdim));
        if self.iteration is 0:
            self._initialise()
        self._import_packing_configuration(fname)
        self.potential = HS_WCA(self.eps, self.sca, self.hs_radii, boxvec=self.boxv)
        hess = self.potential.getHessian(self.coords)
        #print hess
        #hess = self.potential.NumericalHessian(self.coords)
        w, v = np.linalg.eig(hess)
        print w
        self.eigenvalues.extend(w)
        #break down hessian into diagonal elements and  compute their eigenvalues
#        for i in xrange(self.nparticles):
#            i1 = self.bdim*i
#            self.hess_block = hess[i1:i1+self.bdim,i1:i1+self.bdim]
#            w, v = np.linalg.eig(self.hess_block)
#            self.eigenvalues.extend(w)
        self.iteration+=1        
    
    def run(self):
        """run generate packings"""
        for fname in os.listdir(self.packings_dir):
            if 'xyzd' in fname:
                print fname
                self.one_iteration(fname)
        self._histogram_eigenvalues()
    
    def _print_initialise(self):
        base_directory = self.base_directory
        trymakedir(base_directory)
        self._print_parameters()
    
    def _print(self, n):
        """print eigenvalues"""
    
    def _print_parameters(self):
        """writes the simulation parameters"""
        fname = '{}/analyse_jammed_packing.config'.format(self.base_directory)
        f = open(fname,'w')
        f.write('#AUTOMATICALLY GENERATED FILE - DO NOT MODIFY BY HAND\n')
        f.write('#Explore_Jammed_Packings wrapper class input parameters\n')
        f.write('[IMPORTED_JAMMED_PACKING]\n')
        f.write('nparticles: {}\n'.format(self.nparticles))
        f.write('packing_fraction: {}\n'.format(self.imp_packing_frac))
        f.write('boxdim: {}\n'.format(self.bdim))
        f.write('ndim: {}\n'.format(self.ndim))
        f.write('boxv: ')
        for val in self.boxv:
            f.write('{} '.format(val))
        f.write('\n')
        assert(self.sca >0)
        f.write('sca: {}\n'.format(self.sca))
        f.close()
    
if __name__ == "__main__":
    
    analyse = _analyse_jammed_packings()
    start=time.time()
    analyse.run()
    end=time.time()
    print "time elapsed",end-start
    
    
        
                
            
              
                
                
                
