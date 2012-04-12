# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import accept_tests.metropolis as metropolis
import copy


class BasinHopping:
  """A class to run the basin hopping algorithm
  
    coords: 
        The initial set of coordinates.  A one dimensional list or numpy array

    potential: 
        A class implimenting the potential.  The class must have the
        following functions implimented

        energy = potential.getEnergy( coords )
        energy, gradient = potential.getEnergyGradient( coords )

    takeStep: 
        The function which randomly perterbs the system, e.g. random
        dispacement.  It takes the form

        takeStep(coords)

    acceptTests:  ([]) 
        An optional list of functions which return False if a quench should be
        rejected.  The Metropolis test is added to this list by default unless
        the input "nometropolis" is set to False. Each test in the list takes
        the form
 
        accept = test(Eold, Enew, new_coords):

    temperature:  (1.0)
        The temperature used in the metropolis criterion.  If no temperature is
        passed, the default 1.0 is used unless the flag "nometropolis" is set
        to False

    nometropolis: (False)
        Flag to disable the Metropolis accept reject test.
    
    event_after_step:  ([])
        An optional list of functions which act just after each monte carlo
        round.  Each even in the list takes the form

        event(Equench_new, newcoords, acceptstep)
  """
  def __init__(self, coords, potential, takeStep, storage=None, event_after_step=[], \
          acceptTests=[],  \
          temperature=1.0, \
          nometropolis=False \
          ):
    #note: make a local copy of lists of events so that an inputted list is not modified.
    self.coords = coords
    self.storage = storage
    self.potential = potential
    self.takeStep = takeStep
    self.event_after_step = copy.copy(event_after_step)
    self.acceptTests = copy.copy(acceptTests)
    self.temperature = temperature
    self.nometropolis = nometropolis

    if not self.nometropolis:
        self.metrop_test = metropolis.Metropolis(self.temperature)
        self.acceptTests.append( self.metrop_test.acceptReject )

  def run(self, nsteps):
    #########################################################################
    #do initial quench
    #########################################################################
    print "doing monteCarlo, nsteps", nsteps
    potel = self.potential.getEnergy(self.coords)
    print "initial energy", potel

    if(self.storage):
      self.storage(potel, self.coords)
      
    potel, V = self.potential.getEnergyGradient(self.coords)
    print "max gradient", np.max(V), potel
    print "minimizing initial coords"

    newcoords, Equench, self.rms, self.funcalls = self.quench(self.coords)
    Equench_new = Equench
    print "Qu  ", 0, "E=", Equench, "quench_steps= ", self.funcalls, "RMS=", self.rms, "Markov E= ", Equench_new

    coords = newcoords

    if(self.storage):
      self.storage(Equench, coords)
    
    for istep in xrange(nsteps):
        acceptstep, newcoords, Equench_new = self.mcStep(self.coords, Equench)
        print "Qu  ", istep+1, "E=", Equench_new, "quench_steps= ", self.funcalls, "RMS=", self.rms, "Markov E= ", Equench, "accepted=", acceptstep
        for event in self.event_after_step:
            event(Equench_new, newcoords, acceptstep)
        if acceptstep:
            if(self.storage):
              self.storage(Equench_new, newcoords)
        self.coords = newcoords
        Equench = Equench_new
    
  def quench(self, coords):
    import scipy.optimize.lbfgsb
    #newcoords, newE = steepest_descent.steepestDescent(potential.getEnergyGradient, coords, 100)
    newcoords, newE, dictionary = scipy.optimize.fmin_l_bfgs_b(self.potential.getEnergyGradient, coords, iprint=-1, pgtol=1e-3)

    V = dictionary["grad"]
    funcalls = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    if warnflag > 0:
        print "warning: problem with quench: ",
        if warnflag == 1:
            print "too many function evaluations"
        else:
            print dictionary['task']

    rms = V.std()
    #print "quench: Ef=", newE, "steps=", funcalls, "max(V)=", np.max(np.abs(V)), "RMS=", rms
    return newcoords, newE, rms, funcalls 
    
  def mcStep(self, coordsold, Equench_old):
    """take one monte carlo basin hopping step"""
    #########################################################################
    #take step
    #########################################################################
    coords = coordsold.copy() #make  a working copy
    self.takeStep(coords)

    #########################################################################
    #quench
    #########################################################################
    qcoords, Equench, self.rms, self.funcalls = self.quench (coords)

    #########################################################################
    #check whether step is accepted with user defined tests.  If any returns
    #false then reject step.
    #########################################################################
    acceptstep = True
    for test in self.acceptTests:
        if not test(Equench_old, Equench, qcoords):
            acceptstep = False

    #########################################################################
    #return
    #########################################################################
    if acceptstep:
        return acceptstep, qcoords, Equench
    else:
        return acceptstep, coordsold, Equench_old

