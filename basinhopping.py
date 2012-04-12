# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import accept_tests.metropolis as metropolis
import copy
import quench


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

    quenchRoutine:  (quench.quench)
        Optionally pass a non-default quench routine.
        
  """

  def __init__(self, coords, potential, takeStep, storage=None, event_after_step=[], \
          acceptTests=[],  \
          temperature=1.0, \
          nometropolis=False, \
          quenchRoutine = quench.quench \
          ):
    #note: make a local copy of lists of events so that an inputted list is not modified.
    self.coords = copy.copy(coords)
    self.storage = storage
    self.potential = potential
    self.takeStep = takeStep
    self.event_after_step = copy.copy(event_after_step)
    self.acceptTests = copy.copy(acceptTests)
    self.temperature = temperature
    self.nometropolis = nometropolis
    self.quenchRoutine = quenchRoutine

    if not self.nometropolis:
        self.metrop_test = metropolis.Metropolis(self.temperature)
        self.acceptTests.append( self.metrop_test.acceptReject )

    self.stepnum = 0

    #########################################################################
    #store intial structure
    #########################################################################
    potel = self.potential.getEnergy(self.coords)
    print "initial energy", potel
    if(self.storage):
        self.storage(potel, self.coords)
      
    #########################################################################
    #do initial quench
    #########################################################################
    newcoords, Equench, self.rms, self.funcalls = self.quenchRoutine(self.coords, self.potential.getEnergyGradient)
    Equench_new = Equench
    print "Qu  ", self.stepnum, "E=", Equench, "quench_steps= ", self.funcalls, "RMS=", self.rms, "Markov E= ", Equench_new
    self.coords = newcoords

    if(self.storage):
        self.storage(Equench, self.coords)

    self.myE = Equench
    
  def run(self, nsteps):

    Equench = self.myE
    for istep in xrange(nsteps):
        self.stepnum += 1
        acceptstep, newcoords, Equench_new = self.mcStep(self.coords, Equench)
        print "Qu  ", self.stepnum, "E=", Equench_new, "quench_steps= ", self.funcalls, "RMS=", self.rms, "Markov E= ", Equench, "accepted=", acceptstep
        for event in self.event_after_step:
            event(Equench_new, newcoords, acceptstep)
        if acceptstep:
            if(self.storage):
              self.storage(Equench_new, newcoords)
        self.coords = newcoords
        Equench = Equench_new
    
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
    qcoords, Equench, self.rms, self.funcalls = self.quenchRoutine (coords, self.potential.getEnergyGradient)

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

