# -*- coding: iso-8859-1 -*-
import numpy as np
import multiprocessing as mp


def getTemps(Tmin, Tmax, nreplicas):
    """
    set up the temperatures
    distribute them exponentially
    """
    #dT = (Tmax - Tmin) / (nreplicas-1)
    CTE = np.exp( np.log( Tmax / Tmin ) / (nreplicas-1) )
    Tlist = [Tmin* CTE**i for i in range(nreplicas)]
    return Tlist


class MCProcess(mp.Process):
    def __init__(self, mcsys, conn):
        mp.Process.__init__(self)
        self.mcsys = mcsys
        self.conn = conn
    
    def run(self):
        #this redefines mp.Process.run
        #self.mcsys.run(50000)
        while 1:
            message = self.conn.recv()
            #print "message", message
            if message[0] == "kill":
                print "terminating", self.name
                return
            elif message[0] == "return system":
                #sending the system doesn't work
                self.conn.send(self.mcsys)
            elif message[0] == "return energy temperature":
                self.conn.send((self.mcsys.markovE, self.mcsys.temperature))
            elif message[0] == "run":
                nsteps = message[1]
                self.mcsys.run(nsteps) 
                #print self.name, "finished", nsteps, "steps"
                self.conn.send("done")
            elif message[0] == "return energy coords":
                #print self.name, "sending energy coords"
                self.conn.send((self.mcsys.markovE, self.mcsys.coords))
            elif message[0] == "replace energy coords":
                #self.conn.send((self.mcsys.markovE, self.mcsys.coords))
                #print self.name, message[0], self.mcsys.markovE, message[1]
                self.mcsys.markovE = message[1]
                self.mcsys.coords= message[2].copy()
            else:
                print "unknown message:"
                print message


class PTMC(object):
    """A class to run Monte Carlo parallel tempering

    eventually I'd like to submit one mc/basinhopping class and use deepcopy
    to make the replicas.  But it's not working for me.  
    """
    #def __init__(self, mcobject, Tmin = 1., Tmax = 1.2, nreplicas = 4  ):
    def __init__(self, replicas ):
        self.replicas = replicas
        self.nreplicas = len(self.replicas)
        self.exchange_frq = 100
        self.step_num = 0
        
        self.ex_outstream = open("exchanges", "w")
        
        self.replicas_par = []
        self.communicators = []
        for rep in self.replicas:
            parent_conn, child_conn = mp.Pipe()
            rep_par = MCProcess( rep, child_conn)
            self.replicas_par.append( rep_par )
            self.communicators.append( parent_conn )
            rep_par.start()

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
    def getSystems(self):
        #this function doesn't work
        self.replicas_final = []
        for conn in self.communicators:
            conn.send(("return system",))
            rep = conn.recv()
            print "getSystems E", rep.markovE
            self.replicas_final.append(rep)
       
    def end(self):
        #self.getSystems()
        for conn in self.communicators:
            conn.send(("kill",))
        for rep in self.replicas_par:
            rep.join()
            
    def runNoExchanges(self, nsteps):
        for conn in self.communicators:
            conn.send(("run", nsteps))
        #wait till they're all done
        for conn in self.communicators:
            #print "waiting for 'done'"
            message = conn.recv()
            if message != "done":
                print "runNoExchanges> received unexpected message", message
            #else:
                #print "received message", message
         

    def run(self, nsteps):
        stepnum = 0
        while stepnum < nsteps:
            self.runNoExchanges(self.exchange_frq)
            stepnum += self.exchange_frq
            self.step_num += self.exchange_frq
            self.tryExchangePar()
            
    def doExchangePar(self, k):
        """
        do parallel tempering exchange between replicas k and k+1
        
        should we exchange coords or temperature?  Exchanging temperature is faster.  
        Exchanging coords is probably simpler.
        """
        print "exchanging", k, k+1
        self.communicators[k].send(("return energy coords",))
        E1, coords1 = self.communicators[k].recv()
        
        self.communicators[k+1].send(("return energy coords",))
        E2, coords2 = self.communicators[k+1].recv()
        
        self.communicators[k].send(("replace energy coords", E2, coords2))
        self.communicators[k+1].send(("replace energy coords", E1, coords1))


    
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

    def getRepEnergyT(self, k):
        self.communicators[k].send(("return energy temperature",))
        return self.communicators[k].recv()


    def tryExchangePar(self):
        #choose which pair to try and exchange
        k = np.random.random_integers( 0, self.nreplicas - 2)
        #print "trying exchange", k, k+1
        #determine if the exchange will be accepted
        E1, T1 = self.getRepEnergyT(k)
        E2, T2 = self.getRepEnergyT(k+1)
        deltaE = E1 - E2
        deltabeta = 1./T1 - 1./T2
        w = min( 1. , np.exp( deltaE * deltabeta ) )
        rand = np.random.rand()
        if w > rand:
            #accept exchange
            self.ex_outstream.write("accepting exchange %d %d %g %g %g %g %d\n" % (k, k+1, E1, E2, T1, T2, self.step_num) )
            self.doExchangePar(k)
        #else:
            #self.ex_outstream.write("rejecting exchange %d %d %g %g\n" % (k, k+1, w, rand) )
            #print "rejecting exchange ", k, k+1, w, rand

        

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
            self.doExchangePar(k)
        #else:
            #print "rejecting exchange ", k, k+1, w, rand
