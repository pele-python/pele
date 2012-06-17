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
  
      acceptTests:  (None) 
          Acceptance criterion for check.  If None is given, metropolis is used.
  
      confCheck: ([])
          list of checks if current configuration is valid. This is executed before acceptTest
          and accepTest is only called if all checks return True.
  
      temperature:  (1.0)
          The temperature used in the metropolis criterion.  If no temperature is
          passed, the default 1.0 is used unless the flag "nometropolis" is set
          to False

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
            acceptTest=None,  \
            temperature=1.0, \
            quenchRoutine = quench.quench, \
            quenchParameters = dict(), \
            confCheck = [], \
            outstream = sys.stdout
            ):
        #########################################################################
        #initialize MonteCarlo base class
        #########################################################################
        MonteCarlo.__init__(self, coords, potential, takeStep, \
                            storage=storage, \
                            event_after_step=event_after_step, \
                            acceptTest=acceptTest,  \
                            temperature=temperature, \
                            confCheck = confCheck, \
                            outstream=outstream)

        self.quenchRoutine = quenchRoutine
        self.quenchParameters = quenchParameters
        
        #########################################################################
        #do initial quench
        #########################################################################
        self.markovE_old = self.markovE
        newcoords, Equench, self.rms, self.funcalls = \
            self.quenchRoutine(self.coords, self.potential.getEnergyGradient, **self.quenchParameters)

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
        self.takeStep.takeStep(self.coords_after_step, driver=self)

        #########################################################################
        #quench
        #########################################################################
        ret = self.quenchRoutine(self.coords_after_step, \
                                 self.potential.getEnergyGradient, **self.quenchParameters)
        self.trial_coords = ret[0]
        self.trial_energy = ret[1] 
        self.rms = ret[2]
        self.funcalls = ret[3]

        #########################################################################
        # check if step is a valid configuration, otherwise reject
        #########################################################################
        self.acceptstep = True
        for check in self.confCheck:
            if not check(self.trial_energy, self.trial_coords, driver=self):
                self.acceptstep=False
        
        #########################################################################
        #check whether step is accepted with user defined tests.  If any returns
        #false then reject step.
        #########################################################################
        if self.acceptstep:
            self.acceptstep = self.acceptTest(self.markovE, self.trial_energy, self.coords, self.trial_coords)

        #########################################################################
    
        
        
        #return new coords and energy and whether or not they were accepted
        #########################################################################
        return self.acceptstep, self.trial_coords, self.trial_energy


    def printStep(self):
        if self.stepnum % self.printfrq == 0:
            if self.outstream != None:
                self.outstream.write( "Qu   " + str(self.stepnum) + " E= " + str(self.trial_energy) + " quench_steps= " + str(self.funcalls) + " RMS= " + str(self.rms) + " Markov E= " + str(self.markovE_old) + " accepted= " + str(self.acceptstep) + "\n" )
    
    def __getstate__(self):
        ddict = self.__dict__.copy();
        del ddict["outstream"]
        del ddict["potential"]
        del ddict["acceptTests"]
        return ddict #.items()
    
    def __setstate__(self, dct):
        self.__dict__.update(dct)
        self.outstream = sys.stdout
        self.acceptTests = [self.metrop_test.acceptReject]
                
