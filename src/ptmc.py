# -*- coding: iso-8859-1 -*-
import numpy as np


def getTemps(Tmin, Tmax, nreplicas):
    """
    set up the temperatures
    distribute them exponentially
    """
    #dT = (Tmax - Tmin) / (nreplicas-1)
    CTE = np.exp( np.log( Tmax / Tmin ) / (nreplicas-1) )
    Tlist = [Tmin* CTE**i for i in range(nreplicas)]
    return Tlist




class PTMC(object):
    """A class to run basin hopping parallel tempering

    eventually I'd like to submit one mc/basinhopping class and use deepcopy
    to make the replicas.  But it's not working for me.  
    """
    #def __init__(self, mcobject, Tmin = 1., Tmax = 1.2, nreplicas = 4  ):
    def __init__(self, replicas ):
        self.replicas = replicas
        self.nreplicas = len(self.replicas)
        self.exchange_frq = 10
        
        self.ex_outstream = open("exchanges", "w")

        """
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
        """

        #############################################################
        #set up the replicas
        """
        We must be very careful here when we initialize multiple instances.
        They must each be completely independent.  e.g. must make copies
        of the classes called by mcobject
        """
        #############################################################
        """
        self.replicas = []
        for i in range(self.nreplicas):
            T = self.Tlist[i]
            replica = copy.deepcopy(mcobject)
            self.replica.append( replica )
        """
        

    def run(self, nsteps):
        stepnum = 0
        while stepnum < nsteps:
            for rep in self.replicas:
                rep.run( self.exchange_frq )
            stepnum += self.exchange_frq
            self.tryExchange()
            
    def doExchange(self, k):
        """
        do parallel tempering exchange
        
        should we exchange coords or temperature?  Exchanging temperature is faster.  
        Exchanging coords is probably simpler.
        """
        E1 = self.replicas[k].markovE
        coords1 = np.copy( self.replicas[k].coords )
        self.replicas[k].markovE = self.replicas[k+1].markovE 
        self.replicas[k].coords = np.copy( self.replicas[k+1].coords )
        self.replicas[k+1].markovE = E1
        self.replicas[k+1].coords = np.copy(coords1)


    def tryExchange(self):
        #choose which pair to try and exchange
        k = np.random.random_integers( 0, self.nreplicas - 2)
        #print "trying exchange", k, k+1
        #determine if the exchange will be accepted
        deltaE = self.replicas[k].markovE - self.replicas[k+1].markovE
        deltabeta = 1./self.replicas[k].temperature - 1./self.replicas[k+1].temperature
        w = min( 1. , np.exp( deltaE * deltabeta ) )
        rand = np.random.rand()
        if w > rand:
            #accept step
            self.ex_outstream.write("accepting exchange %d %d %g %g\n" % (k, k+1, w, rand) )
            self.doExchange(k)
        #else:
            #print "rejecting exchange ", k, k+1, w, rand
