# -*- coding: iso-8859-1 -*-
from __future__ import absolute_import
import sys
from .accept_tests import metropolis as metropolis
import copy
import numpy as np
from pele.optimize import Result


class MonteCarlo(object):
    """A class to run the Monte Carlo algorithm

    Parameters
    ----------
    coords : numpy array, one dimensional
        The initial set of coordinates.
    potential : potential object
        A class implementing the potential.  The class must have the
        following functions implemented::
    
            potential.getEnergy(coords)
            potential.getEnergyGradient(coords)
    
    takeStep : object
        The class which randomly perterbs the system, e.g. random
        displacement.  It must have two functions implemented which
        are called like::
    
            takeStep.takeStep(coords, driver=self)       #  takes the step
            takeStep.updateStep(acceptstep, driver=self) #  for adaptive step size management
    
    acceptTest : callable, optional
        Acceptance criterion for monte carlo move.  If None is given, metropolis is used.
        It must have the form::

            bool = acceptTest(Eold, Enew, old_coords, new_coords)
        
    confCheck : list of callables, optional
        list of checks if current configuration is valid. This is executed before acceptTest
        and accepTest is only called if all checks return True.  The checks are called like::
        
            check(trial_energy, trial_coords, driver=self)
        
        and must return a bool
    
    temperature : float, optional
        The temperature used in the metropolis criterion.  
    event_after_step :  list of callables, optional
        these are called just after each monte carlo
        round.  Each event in the list takes the form:::
    
            event(energy, coords, acceptstep)
    
    outstream : open file object, optional
        The file stream to print quench information to.  None for no printing.
        Default to standard out.
    store_initial : bool, optional
        if True store initial structure
    
    See Also
    --------
    pele.potentials, pele.takestep, pele.storage, pele.accept_tests
    """
    
    insert_rejected = False
  
    def __init__(self, coords, potential, takeStep, storage=None, event_after_step=None, acceptTest=None,
                 temperature=1.0, confCheck=None, outstream=sys.stdout, store_initial=True, iprint=1):
        # note: make a local copy of lists of events so that an inputted list is not modified.
        if confCheck is None: confCheck = []
        if event_after_step is None: event_after_step = []
        self.coords = np.copy(coords)
        self.storage = storage
        self.potential = potential
        self.takeStep = takeStep
        self.event_after_step = copy.copy(event_after_step) # not deepcopy
        self.temperature = temperature
        self.naccepted = 0

        self.result = Result()
        self.result.nfev = 0
        
        self.outstream = outstream
        self.printfrq = iprint # controls how often printing is done
        self.confCheck = confCheck
    
        if acceptTest:
            self.acceptTest = acceptTest 
        else:
            self.acceptTest = metropolis.Metropolis(self.temperature)
        
        self.stepnum = 0
    
        #########################################################################
        # store intial structure
        #########################################################################
        energy = self.potential.getEnergy(self.coords)
        self.result.nfev += 1
        if self.storage and store_initial:
            self.storage(energy, self.coords)
          
        self.markovE = energy
        
        self.result.energy = self.markovE
        self.result.coords = self.coords.copy()

    
    def setPrinting(self, ostream="default", frq=None):
        """change how the printing is done
        
        Parameters
        -----------
        ostream : open file or None
            where to print data
        frq : int or None
            how often to print data.  if None, don't change current value
            
        """
        if ostream != "default": # ostream = None is valid input
            self.outstream = ostream
        if frq is not None:
            self.printfrq = frq
    
    def addEventAfterStep(self, event):
        """add an even to the list event_after_step """
        self.event_after_step.append( event )
            
    
    def _mcStep(self):
        """take one monte carlo basin hopping step
        """
        self.trial_coords = self.coords.copy() # make  a working copy
        #########################################################################
        # take step
        #########################################################################
        self.takeStep.takeStep(self.trial_coords, driver=self)
                
        #########################################################################
        # calculate new energy
        #########################################################################
        self.trial_energy = self.potential.getEnergy(self.trial_coords)
        self.result.nfev += 1
        
        
        
        #########################################################################
        # check if step is a valid configuration, otherwise reject
        #########################################################################
        self.acceptstep = True
        self.config_ok = True
        for check in self.confCheck:
            if not check(self.trial_energy, self.trial_coords, driver=self):
                self.acceptstep=False
                self.config_ok = False
        
        #########################################################################
        # check whether step is accepted with user defined tests.  If any returns
        # false then reject step.
        #########################################################################
        if self.acceptstep:
            self.acceptstep = self.acceptTest(self.markovE, self.trial_energy,
                                              self.coords, self.trial_coords)
            
        #########################################################################
        # return new coords and energy and whether or not they were accepted
        #########################################################################
        return self.acceptstep, self.trial_coords, self.trial_energy
  
    def run(self, nsteps):
        """do multiple iterations"""
        # take nsteps
        for istep in range(nsteps):
            self.takeOneStep()
    
    def takeOneStep(self):
        """one cycle of the routine
        """
        self.stepnum += 1
        self.markovE_old = self.markovE
        acceptstep, newcoords, newE = self._mcStep()
        self.printStep()
        if self.storage and (self.insert_rejected or acceptstep) and self.config_ok:
            self.storage(newE, newcoords)

        if acceptstep:
            self.coords = newcoords
            self.markovE = newE
            self.naccepted += 1
            
            if self.markovE < self.result.energy:
                self.result.energy = self.markovE
                self.result.coords = self.coords.copy()

        self.takeStep.updateStep(acceptstep, driver=self)
        for event in self.event_after_step:
            event(self.markovE, self.coords, acceptstep)

    def printStep(self):
        if self.stepnum % self.printfrq == 0:
            if self.outstream != None:
                self.outstream.write("MCstep    %12d  E= %20.12g  markov E= %20.12g accepted= %s\n" 
                        % (self.stepnum, self.trial_energy, self.markovE_old, str(self.acceptstep) )  )



if __name__ == "__main__":
    from pele.systems import LJCluster
    natoms = 13
    sys = LJCluster(natoms)
    pot = sys.get_potential()
    coords = sys.get_random_configuration()
    takestep = sys.get_takestep(stepsize=.1)
    mc = MonteCarlo(coords, pot, takestep)
    mc.run(100)

