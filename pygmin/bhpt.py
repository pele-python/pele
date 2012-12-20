# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import accept_tests.metropolis as metropolis
import copy
from optimize import mylbfgs
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

      quenchRoutine:  (mylbfgs)
          Optionally pass a non-default quench routine.

    """
    def __init__(self, coords, potential, takeStep, storage=None, event_after_step=[], \
            acceptTests=[],  \
            nometropolis=False, \
            quenchRoutine = mylbfgs, \
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

        self.nreplicas = nreplicas
        self.exchange_frq = 10

        #set up the temperatures
        #distribute them exponentially
        dT = (Tmax - Tmin) / (self.nreplicas-1)
        CTE = np.exp( np.log( Tmax / Tmin ) / (self.nreplicas-1) )
        self.Tlist = [Tmin* CTE**i for i in range(self.nreplicas)]
        print "Tlist", self.Tlist

        self.streams = []
        #set up the outstreams
        for i in range(self.nreplicas):
            self.streams.append( open("minGMIN_out." + str(i), "w" ) )

        #############################################################
        #set up the replicas
        """
        We must be very careful here when we initialize multiple instances of
        BasinHopping by passing classes and routines (e.g. takeStep).  For some
        classes we can pass a reference to the same class instance.  For others
        it is very important that each BasinHopping instance has it's own copy
        of the class.
        """
        #############################################################
        self.replicas = []
        for i in range(self.nreplicas):
            T = self.Tlist[i]
            replica = bh.BasinHopping( self.coords, self.potential, self.takeStep, temperature = T, outstream = self.streams[i])
            self.replicas.append( replica )

    def run(self, nsteps):

        for istep in xrange(nsteps/self.exchange_frq):
            for rep in self.replicas:
                rep.run( self.exchange_frq )
            self.tryExchange()

    def tryExchange(self):
        k = np.random.random_integers( 0, self.nreplicas - 2)
        print "trying exchange", k, k+1
        deltaE = self.replicas[k].markovE - self.replicas[k+1].markovE
        deltabeta = 1./self.replicas[k].temperature - 1./self.replicas[k+1].temperature
        w = min( 1. , np.exp( deltaE * deltabeta ) )
        rand = np.random.rand()
        if w > rand:
            #accept step
            print "accepting exchange ", k, k+1, w, rand
            E1 = self.replicas[k].markovE
            coords1 = copy.copy( self.replicas[k].coords )
            self.replicas[k].markovE = self.replicas[k+1].markovE 
            self.replicas[k].coords = copy.copy( self.replicas[k+1].coords )
            self.replicas[k+1].markovE = E1
            self.replicas[k+1].coords = coords1
        else:
            print "rejecting exchange ", k, k+1, w, rand
