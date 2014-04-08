from __future__ import division
import numpy as np
import abc
import os
from pele.potentials import HS_WCA
from playground.basin_volume import BV_MCrunner
from utils import *
import ConfigParser
import time

class _configure_bv_mcrunner(object):
    """
    this is an abstract class that implements the basic components of a configure bv_mcrunner class,
    and declares a number of abstract methods which should be implemented in all inheriting classes
    *nparticles: number of particles
    *bdim: dimensionality of the box
    *ndim: dimensionality of the problem (i.e. size of the coordinates array)
    *packing_frac: target jammed packing fraction
    *boxv: an array of size bdim that contains the vectors defining the box
    """
        
    def __call__(self, fname, k=100.0, temperature=1.0, niter=1e4, stepsize=1e-2, dtol=1e-4, eps=1., hmin=0, 
                 hmax=3, hbinsize=1e-2, acceptance=0.2, adjustf=0.9, adjustf_niter = 2000, adjustf_navg = 100, 
                 opt_dtmax=1, opt_maxstep=0.5, opt_tol=1e-3, opt_nsteps=1e5, packings_dir='jammed_packings'):
        dname = fname
        if dname.endswith('.xyzd'):
            dname = dname[:-5]
        self.base_directory = os.path.join(os.getcwd(),'explore_bv_'+str(dname))
        self.packings_dir = os.path.join(os.getcwd(),packings_dir)
        self.configpath = os.path.join(packings_dir,'jammed_packings.config')
        self.fname = fname
        #self.mc_params = dict(k=k, temperature=temperature, )
        self.mc_params = {'k':k,'temperature':temperature,'niter':niter,'stepsize':stepsize,'dtol':dtol,'eps':eps,'hmin':hmin,'hmax':hmax,
                      'hbinsize':hbinsize,'acceptance':acceptance,'adjustf':adjustf,'adjustf_niter':adjustf_niter,'adjustf_navg':adjustf_navg,
                      'opt_dtmax':opt_dtmax,'opt_maxstep':opt_maxstep,'opt_tol':opt_tol,'opt_nsteps':opt_nsteps}
        self._initialise()
        #construct mcrunner
        #self.coords is origin, set initial configuration and origin to be the same
        mcrunner = BV_MCrunner(self.coords, self.coords, self.hs_radii, self.boxv, self.sca, rattlers=self.rattlers, k=k, 
                               temperature=temperature, niter=niter, stepsize=stepsize, dtol=dtol, eps=eps, hmin=hmin, 
                               hmax=hmax, hbinsize=hbinsize, acceptance=acceptance, adjustf=adjustf, adjustf_niter = adjustf_niter, 
                               adjustf_navg = adjustf_navg, opt_dtmax=opt_dtmax, opt_maxstep=opt_maxstep, opt_tol=opt_tol,
                               opt_nsteps=opt_nsteps)
        return mcrunner 
        
    def _initialise(self):
        """initialisation function"""
        self._import_packing_config_file()
        self._import_packing_configuration()
        #here there should be the import rattlers function
        self.rattlers = None
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
        
    def _import_packing_configuration(self):
        """imports the coordinates and data relative to the shinitape of the particles
            this should be run in initialise()
        """
        path = os.path.join(self.packings_dir,self.fname)
        self.coords, self.hs_radii = read_xyzd(path)
            
    def _print_initialise(self):
        base_directory = self.base_directory
        trymakedir(base_directory)
        self._print_parameters()
    
    def _print_parameters(self):
        """writes the simulation parameters"""
        dname = self.fname 
        if dname.endswith('.xyzd'):
            dname = dname[:-5]
        fname = '{}/explore_{}.config'.format(self.base_directory,dname)
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
        f.write('[MCRUNNER]\n')
        for key, value in self.mc_params.iteritems() :
            f.write('{}: {}\n'.format(key,value)) 
        f.close()
    
if __name__ == "__main__":
    
    sim = _configure_bv_mcrunner()
    mcrunner = sim('jammed_packing0.xyzd')
    start=time.time()
    mcrunner.run()
    end=time.time()
    print end-start
    status = mcrunner.get_status()
    print status
    mcrunner.show_histogram()
    
    
        
                
            
              
                
                
                
