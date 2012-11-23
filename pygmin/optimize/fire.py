'''
Created on 30 Apr 2012

@author: ruehle
'''

import numpy as np
import math


class Fire(object):
    '''
    classdocs
    '''
    
    def __init__(self, coords, potential, restart=None, logfile='-', trajectory=None,
                 dt=0.1, maxmove=0.5, dtmax=1., Nmin=5, finc=1.1, fdec=0.5,
                 astart=0.1, fa=0.99, a=0.1, iprint=-1,
                 alternate_stop_criterion = None):
        #Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        self.dt = dt
        self.Nsteps = 0
        self.maxmove = maxmove
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
        self.coords = coords
        self.potential=potential
        self.v = None
        self.nsteps=0
        self.iprint = iprint
        self.alternate_stop_criterion = alternate_stop_criterion
        
    def initialize(self):
        self.v = None
        
    def step(self,f):
        coords = self.coords
        if self.v is None:
            self.v = np.zeros((len(coords)))
        else:
            vf = np.vdot(f, self.v)
            if vf > 0.0:
                self.v = (1.0 - self.a) * self.v + self.a * f / np.sqrt(
                    np.vdot(f, f)) * np.sqrt(np.vdot(self.v, self.v))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt * self.finc, self.dtmax)
                    self.a *= self.fa
                self.Nsteps += 1
            else:
                self.v[:] *= 0.0
                self.a = self.astart
                self.dt *= self.fdec
                self.Nsteps = 0

        self.v += self.dt * f
        dr = self.dt * self.v
        if False: #how do we determine maxmove?
            normdr = np.sqrt(np.vdot(dr, dr))
        else:
            normdr = max(np.abs(dr))
        #print "aa",normdr
        if normdr > self.maxmove:
            dr = self.maxmove * dr / normdr
        self.coords= coords + dr
     
    def run(self, fmax=1e-3, steps=100000):
        """Run structure optimization algorithm.

        This method will return when the forces on all individual
        atoms are less than *fmax* or when the number of steps exceeds
        *steps*."""

        self.fmax = fmax
        step = 0
        while step < steps:
            E,f = self.potential(self.coords)
            #self.call_observers()
            #print E
            if self.alternate_stop_criterion is None:
                i_am_done = self.converged(f)
            else:
                i_am_done = self.alternate_stop_criterion(energy = E, gradient = f, 
                                                          tol=self.fmax)                
            if i_am_done:
                return
            self.step(-f)
            self.nsteps += 1
            if self.iprint > 0:
                if step % self.iprint == 0:
                    rms = np.linalg.norm(f)/np.sqrt(len(f))
                    print "fire:", step, E, rms
            step += 1
                    

    def converged(self, forces=None):
        """Did the optimization converge?"""
        #if forces is None:
        #    forces = self.atoms.get_forces()
        #if hasattr(self.atoms, 'get_curvature'):
        #    return (forces**2).sum(axis=1).max() < self.fmax**2 and \
        #           self.atoms.get_curvature() < 0.0
        #print forces
        #print (forces**2).sum()
        #print max(forces)
        #return (forces**2).sum().max() < self.fmax**2
        #print np.linalg.norm(forces)/math.sqrt(len(forces))
        return np.linalg.norm(forces)/math.sqrt(len(forces)) < self.fmax

if __name__ == "__main__":
    import pygmin.potentials.lj as lj
    pot = lj.LJ()
    coords = 10.*np.random.random(300)
    opt = Fire(coords, pot.getEnergyGradient, dtmax=0.1, dt=0.01, maxmove=0.01)
    opt.run(fmax=1e-1,steps=10000)
    print(pot.getEnergy(opt.coords))
    print opt.nsteps
