# -*- coding: iso-8859-1 -*-
import sys
import optimize.quench as quench
from mc import MonteCarlo

class BasinHopping(MonteCarlo):
    """A class to run the basin hopping algorithm

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

      quenchRoutine:  (quench.quench)
          Optionally pass a non-default quench routine.

      outstream: (stdout)
          the file stream to print quench information to
    """

    def __init__(self, coords, potential, takeStep, storage=None, event_after_step=[], \
            acceptTests=[],  \
            temperature=1.0, \
            nometropolis=False, \
            quenchRoutine = quench.quench, \
            outstream = sys.stdout
            ):
        #########################################################################
        #initialize MonteCarlo base class
        #########################################################################
        MonteCarlo.__init__(self, coords, potential, takeStep, storage, event_after_step, \
              acceptTests,  \
              temperature, \
              nometropolis, \
              outstream)

        self.quenchRoutine = quenchRoutine

        #########################################################################
        #do initial quench
        #########################################################################
        self.markovE_old = self.markovE
        newcoords, Equench, self.rms, self.funcalls = self.quenchRoutine(self.coords, self.potential.getEnergyGradient)

        self.coords = newcoords
        self.markovE = Equench

        if(self.storage):
            self.storage(self.markovE, self.coords)
        
        #print the initial quench
        self.acceptstep = True
        self.trial_energy = self.markovE
        self.printStep()


    def mcStep(self):
        """
        take one monte carlo basin hopping step

        redefine the MonteCarlo base class step
        """
        self.coords_after_step = self.coords.copy() #make  a working copy
        #########################################################################
        #take step
        #########################################################################
        self.takeStep.takeStep(self.coords_after_step)

        #########################################################################
        #quench
        #########################################################################
        ret = self.quenchRoutine(self.coords_after_step, self.potential.getEnergyGradient)
        self.trial_coords = ret[0]
        self.trial_energy = ret[1] 
        self.rms = ret[2]
        self.funcalls = ret[3]

        #########################################################################
        #check whether step is accepted with user defined tests.  If any returns
        #false then reject step.
        #########################################################################
        self.acceptstep = True
        for test in self.acceptTests:
            if not test(self.markovE, self.trial_energy, self.trial_coords, self.coords_after_step):
                self.acceptstep = False

        #########################################################################
        #return new coords and energy and whether or not they were accepted
        #########################################################################
        return self.acceptstep, self.trial_coords, self.trial_energy


    def printStep(self):
        if self.outstream != None:
            self.outstream.write( "Qu   " + str(self.stepnum) + " E= " + str(self.trial_energy) + " quench_steps= " + str(self.funcalls) + " RMS= " + str(self.rms) + " Markov E= " + str(self.markovE_old) + " accepted= " + str(self.acceptstep) + "\n" )

