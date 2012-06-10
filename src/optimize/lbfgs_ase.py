# -*- coding: utf-8 -*-
import sys
import numpy as np
from optimize_ase import Optimizer
#from ase.utils.linesearch import LineSearch



class Atoms:
    def __init__(self, coords, getEnergy, getGradient, ndim = 3):
        self.coords = coords.copy()
        self.getEnergy = getEnergy
        self.getGradient = getGradient
        self.ndim = ndim
        self.natoms = len(coords) / self.ndim
        self.myshape = [self.natoms, self.ndim]
    def get_potential_energy(self):
        return self.getEnergy(self.coords)
    def get_forces(self):
        return np.reshape( -self.getGradient(self.coords), self.myshape)
    def get_positions(self):
        x = np.reshape(self.coords, self.myshape ).copy()
        #print np.shape(x), self.coords[0:3], x[0,:]
        return x
    def set_positions(self, newcoords):
        x = newcoords.copy()
        x = np.reshape( x, self.natoms*self.ndim )
        self.coords = x.copy()
        #print np.shape(newcoords), self.coords[-3:], newcoords[-1,:], "set"




class LBFGS(Optimizer):
    """Limited memory BFGS optimizer.
    
    A limited memory version of the bfgs algorithm. Unlike the bfgs algorithm
    used in bfgs.py, the inverse of Hessian matrix is updated.  The inverse
    Hessian is represented only as a diagonal matrix to save memory

    """
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, memory=100, damping = 1.0, alpha = 10.0,
                 use_line_search=False):
        """
        Parameters:

        restart: string
            Pickle file used to store vectors for updating the inverse of Hessian
            matrix. If set, file with such a name will be searched and information
            stored will be used, if the file exists.

        logfile: string
            Where should output go. None for no output, '-' for stdout.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        maxstep: float
            How far is a single atom allowed to move. This is useful for DFT
            calculations where wavefunctions can be reused if steps are small.
            Default is 0.04 Angstrom.

        memory: int
            Number of steps to be stored. Default value is 100. Three numpy
            arrays of this length containing floats are stored.

        damping: float
            The calculated step is multiplied with this number before added to
            the positions. 

        alpha: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.
            
        """
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if maxstep is not None:
            if maxstep > 1.0:
                raise ValueError('You are using a much too large value for ' +
                                 'the maximum step size: %.1f Angstrom' % maxstep)
            self.maxstep = maxstep
        else:
            self.maxstep = 0.04

        self.memory = memory
        self.H0 = 1. / alpha  # Initial approximation of inverse Hessian
                            # 1./70. is to emulate the behaviour of BFGS
                            # Note that this is never changed!
        self.damping = damping
        self.use_line_search = use_line_search
        self.p = None
        self.function_calls = 0
        self.force_calls = 0
        print "ase> initial energy", self.func(self.atoms.get_positions())

    def initialize(self):
        """Initalize everything so no checks have to be done in step"""
        self.iteration = 0
        self.s = []
        self.y = []
        self.rho = [] # Store also rho, to avoid calculationg the dot product
                      # again and again

        self.r0 = None
        self.f0 = None
        self.e0 = None
        self.task = 'START'
        self.load_restart = False

    def read(self):
        """Load saved arrays to reconstruct the Hessian"""
        self.iteration, self.s, self.y, self.rho, \
        self.r0, self.f0, self.e0, self.task = self.load()
        self.load_restart = True

    def step(self, f):
        """Take a single step
        
        Use the given forces, update the history and calculate the next step --
        then take it"""
        r = self.atoms.get_positions()
        p0 = self.p
    
        self.update(r, f, self.r0, self.f0)
        
        s = self.s
        y = self.y
        rho = self.rho
        H0 = self.H0

        loopmax = np.min([self.memory, self.iteration])
        a = np.empty((loopmax,), dtype=np.float64)

        ### The algorithm itself:
        # q(3*natoms): forces
        # s[](3*natoms):  list of past changes in positions
        # rho[]:  list of numbers
        # a[]:  list of numbers
        # y[](3*natoms):  list of past changes in gradient
        # H0:  float: initial guess for hessian
        # z(3*natoms):  
        # b: float
        # 
        ### The algorithm itself:
        #see http://en.wikipedia.org/wiki/Limited-memory_BFGS
        q = - f.reshape(-1) 
        for i in range(loopmax - 1, -1, -1):
            a[i] = rho[i] * np.dot(s[i], q)
            q -= a[i] * y[i]
        z = H0 * q
        
        for i in range(loopmax):
            b = rho[i] * np.dot(y[i], z)
            z += s[i] * (a[i] - b)

        self.p = - z.reshape((-1, 3))
        ###
        
        if False:
            print "q", np.shape(q)
            print "s", len(s), np.shape(s)
            print "rho", np.shape(rho)
            print "a", np.shape(a)
            print "a", np.shape(a)
            print "y", np.shape(y)
            print "z", np.shape(z)
            print "self.p", np.shape(self.p)

 
        g = -f
        if self.use_line_search == True:
            e = self.func(r)
            self.line_search(r, g, e)
            dr = (self.alpha_k * self.p).reshape(len(self.atoms),-1)
        else:
            self.force_calls += 1
            self.function_calls += 1
            dr = self.determine_step(self.p) * self.damping
        self.atoms.set_positions(r+dr)
        
        self.iteration += 1
        self.r0 = r
        self.f0 = -g
        self.dump((self.iteration, self.s, self.y, 
                   self.rho, self.r0, self.f0, self.e0, self.task))

    def determine_step(self, dr):
        """Determine step to take according to maxstep
        
        Normalize all steps as the largest step. This way
        we still move along the eigendirection.
        """
        steplengths = (dr**2).sum(1)**0.5
        longest_step = np.max(steplengths)
        if longest_step >= self.maxstep:
            dr *= self.maxstep / longest_step
        
        return dr

    def update(self, r, f, r0, f0):
        """Update everything that is kept in memory

        This function is mostly here to allow for replay_trajectory.
        """
        if self.iteration > 0:
            s0 = r.reshape(-1) - r0.reshape(-1)
            self.s.append(s0)

            # We use the gradient which is minus the force!
            y0 = f0.reshape(-1) - f.reshape(-1)
            self.y.append(y0)
            
            rho0 = 1.0 / np.dot(y0, s0)
            self.rho.append(rho0)

        if self.iteration > self.memory:
            self.s.pop(0)
            self.y.pop(0)
            self.rho.pop(0)


    def replay_trajectory(self, traj):
        """Initialize history from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import PickleTrajectory
            traj = PickleTrajectory(traj, 'r')
        r0 = None
        f0 = None
        # The last element is not added, as we get that for free when taking
        # the first qn-step after the replay
        for i in range(0, len(traj) - 1):
            r = traj[i].get_positions()
            f = traj[i].get_forces()
            self.update(r, f, r0, f0)
            r0 = r.copy()
            f0 = f.copy()
            self.iteration += 1
        self.r0 = r0
        self.f0 = f0

    def func(self, x):
        """Objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.function_calls += 1
        return self.atoms.get_potential_energy()

    def fprime(self, x):
        """Gradient of the objective function for use of the optimizers"""
        self.atoms.set_positions(x.reshape(-1, 3))
        self.force_calls += 1
        # Remember that forces are minus the gradient!
        return - self.atoms.get_forces().reshape(-1)

    def line_search(self, r, g, e):
        print "line search not implemented"
        self.p = self.p.ravel()
        p_size = np.sqrt((self.p **2).sum())
        if p_size <= np.sqrt(len(self.atoms) * 1e-10):
            self.p /= (p_size / np.sqrt(len(self.atoms)*1e-10))
        g = g.ravel()
        r = r.ravel()
        ls = LineSearch()
        self.alpha_k, e, self.e0, self.no_update = \
           ls._line_search(self.func, self.fprime, r, self.p, g, e, self.e0,
                           maxstep=self.maxstep, c1=.23,
                           c2=.46, stpmax=50.)
        if self.alpha_k is None:
            raise RuntimeError("LineSearch failed!")


