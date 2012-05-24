# -*- coding: iso-8859-1 -*-
import numpy as np
from math import *
import sys
import accept_tests.metropolis as metropolis
import copy


class MonteCarlo(object):
    """A class to run the monte carlo algorithm
    
      coords: 
          The initial set of coordinates.  A one dimensional list or numpy array
  
      potential: 
          A class implimenting the potential.  The class must have the
          following functions implimented
  
          energy = potential.getEnergy( coords )
          energy, gradient = potential.getEnergyGradient( coords )
  
      takeStep: 
          The class which randomly perterbs the system, e.g. random
          dispacement.  It must have two functions implemented
  
          takeStep.takeStep(coords)       #actually takes the step
          takeStep.updateStep(coords)     #for adaptive step size management
  
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
  
      outstream: (stdout)
          the file stream to print quench information to
    """
  
    def __init__(self, coords, potential, takeStep, storage=None, event_after_step=[], \
            acceptTests=[],  \
            temperature=1.0, \
            nometropolis=False, \
            outstream = sys.stdout
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
        self.naccepted = 0
    
        self.outstream = outstream
    
        if not self.nometropolis:
            self.metrop_test = metropolis.Metropolis(self.temperature)
            self.acceptTests.append( self.metrop_test.acceptReject )
    
        self.stepnum = 0
    
        #########################################################################
        #store intial structure
        #########################################################################
        energy = self.potential.getEnergy(self.coords)
        if(self.storage):
            self.storage(energy, self.coords)
          
        self.markovE = energy
        
    def mcStep(self):
        """take one monte carlo basin hopping step"""
        self.trial_coords = self.coords.copy() #make  a working copy
        #########################################################################
        #take step
        #########################################################################
        self.takeStep.takeStep(self.trial_coords)
        
        #########################################################################
        #calculate new energy
        #########################################################################
        self.trial_energy = self.potential.getEnergy(self.trial_coords)
        
        
        #########################################################################
        #check whether step is accepted with user defined tests.  If any returns
        #false then reject step.
        #########################################################################
        self.acceptstep = True
        for test in self.acceptTests:
            if not test(self.markovE, self.trial_energy, self.coords, self.trial_coords):
                self.acceptstep = False
    
        #########################################################################
        #return new coords and energy and whether or not they were accepted
        #########################################################################
        return self.acceptstep, self.trial_coords, self.trial_energy
  
    def run(self, nsteps):
        #take nsteps
        for istep in xrange(nsteps):
            self.takeOneStep()
    
    def takeOneStep(self):
        self.stepnum += 1
        self.markovE_old = self.markovE
        acceptstep, newcoords, newE = self.mcStep()
        #self.outstream.write( "Qu   " + str(self.stepnum) + " E= " + str(newE) + " quench_steps= " + str(self.funcalls) + " RMS= " + str(self.rms) + " Markov E= " + str(self.markovE) + " accepted= " + str(acceptstep) + "\n" )
        self.printStep()
        try:
            self.takeStep.updateStep(acceptstep)
        except:
            print "WARNING: takeStep.updateStep() not implemented"
        if acceptstep:
            if(self.storage):
                self.storage(newE, newcoords)
            self.coords = newcoords
            self.markovE = newE
            self.naccepted += 1
        for event in self.event_after_step:
            event(self.markovE, self.coords, acceptstep)

    def printStep(self):
        if self.outstream != None:
            self.outstream.write( "MCstep    %12d  E= %20.12g  markov E= %20.12g accepted= %s\n" % (self.stepnum, self.trial_energy, self.markovE_old, str(self.acceptstep) )  )




