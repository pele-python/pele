# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import accept_tests.metropolis as metropolis
import copy
import quench
import basinhopping as bh


class BHPT:
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
            nometropolis=False, \
            quenchRoutine = quench.quench, \
            Tmin = 1., Tmax = 1.2, nreplicas = 4 \
            ):
        #note: make a local copy of lists of events so that an inputted list is not modified.
        self.coords = coords
        self.storage = storage
        self.potential = potential
        self.takeStep = takeStep
        self.event_after_step = copy.copy(event_after_step)
        self.acceptTests = copy.copy(acceptTests)
        self.nometropolis = nometropolis
        self.quenchRoutine = quenchRoutine

        dT = (Tmax - Tmin) / (nreplicas-1)
        self.Tlist = [Tmin + i*dT for i in range(nreplicas)]
        print "Tlist", self.Tlist

        self.replicas = []
        for T in self.Tlist:
            replica = bh.BasinHopping( self.coords, self.potential, self.takeStep, temperature = T)
            self.replicas.append( replica )

    def run(self, nsteps):

        for istep in xrange(nsteps):
            for rep in self.replicas:
                rep.run(1)

