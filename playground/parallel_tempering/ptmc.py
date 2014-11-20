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
    """
    ****This is still in testing, and definitely not ready for any production runs****
    
    A class to run Monte Carlo parallel tempering
    """
    #def __init__(self, mcobject, Tmin = 1., Tmax = 1.2, nreplicas = 4  ):
    def __init__(self, replicas ):
        self.replicas = replicas
        self.nreplicas = len(self.replicas)
        self.exchange_frq = 100
        self.step_num = 0
        self.use_independent_exchange = False
        
        self.ex_outstream = open("exchanges", "w")
        
        self.replicas_par = []
        self.communicators = []
        for rep in self.replicas:
            parent_conn, child_conn = mp.Pipe()
            rep_par = MCProcess( rep, child_conn)
            self.replicas_par.append( rep_par )
            self.communicators.append( parent_conn )
            rep_par.start()

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
            self.tryExchange()
            
    def doExchangePar(self, k1, k2):
        """
        do parallel tempering exchange between replicas k and k+1
        
        should we exchange coords or temperature?  Exchanging temperature is faster.  
        Exchanging coords is probably simpler.
        """
        #print "exchanging", k1, k2
        self.communicators[k1].send(("return energy coords",))
        E1, coords1 = self.communicators[k1].recv()
        
        self.communicators[k2].send(("return energy coords",))
        E2, coords2 = self.communicators[k2].recv()
        
        self.communicators[k1].send(("replace energy coords", E2, coords2))
        self.communicators[k2].send(("replace energy coords", E1, coords1))


    
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
            self.doExchangePar(k, k+1)
        #else:
            #self.ex_outstream.write("rejecting exchange %d %d %g %g\n" % (k, k+1, w, rand) )
            #print "rejecting exchange ", k, k+1, w, rand

    def tryExchange(self):
        if self.use_independent_exchange:
            self.tryExchangeIndependent()
        else:
            self.tryExchangePar()

    def tryExchangeNoParallel(self):
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
            self.doExchangePar(k, k+1)
        #else:
            #print "rejecting exchange ", k, k+1, w, rand

    def doMultipleExchanges(self, newindices):
        newindices = [ val for val in enumerate(newindices) if val[0] != val[1] ]
        energy_coords = []#[ None for i in range(len(newindices))]
        #get old coordinates and energies
        for newindex, oldindex in newindices:
            self.communicators[oldindex].send(("return energy coords",))
            E, coords = self.communicators[oldindex].recv()
            energy_coords.append( (newindex, E, coords) )
        #replace the energies
        for newindex, E, coords in energy_coords:
            self.communicators[newindex].send(("replace energy coords", E, coords))

             

            

    def tryExchangeIndependent(self):
        ET = [ self.getRepEnergyT(k) for k in range(self.nreplicas) ]
        energies = [ val[0] for val in ET ]
        beta = [ 1./val[1] for val in ET ]
        ptexch = PTExchangeIndependent(beta, energies)
        ptexch.run()
        self.doMultipleExchanges(ptexch.replicas)
        if True:
            ostr = "%d multiple exchanges" % (self.step_num)
            for i in ptexch.replicas:
                ostr += " %d" % (i)
            ostr += "\n"
            self.ex_outstream.write(ostr)
#        for i1, i2 in enumerate(ptexch.replicas):
#            if i1 != i2:
#                self.ex_outstream.write("accepting exchange %d %d %g %g %d\n" % (i1, i2, energies[i1], energies[i2], self.step_num) )
#                self.doExchangePar(i1, i2)


class PTExchangeIndependent(object):
    """
    implement the independence sampling parallel tempering exchange step
    
    this is a replacement for simple neighbor exchange.  
    
    the basic idea is rather than do just one exchange, do many.  Enough so that
    the replicas are distributed independently given the appropriate probability 
    distribution, with no memory of their original positions.
    
    see J. D. Chodera, and M. R. Shirts, JCP 2011
    http://dx.doi.org/10.1063/1.3660669
    """
    def __init__(self, beta, energies):
        self.beta  = np.array(beta)
        self.energies = np.array(energies)
        self.nreps = len(self.beta)
        self.replicas = np.array(range(self.nreps))
        self.nsteps = self.nreps**3 #this is rule of thumb
        self.naccept = 0
        self.acccount = np.zeros(self.nreps)
        """
        note:
        replicas[i] is the label of the replica which is currently at i.  E.g.
        after the monte carlo run.  Replica replicas[i] should be moved to
        position i

        """

    def run(self):
        nsteps = self.nreps**3
        for i in range(nsteps):
            self.mcstep()
        #print self.replicas
        print "accept ratio", float(self.naccept) / nsteps, nsteps, self.naccept
        print self.acccount

    def getPair(self):
        i=0
        j=i
        while i == j:
            i,j = np.random.random_integers(0,self.nreps-1,2)
        return i,j

    def mcstep(self):
        i1, i2 = self.getPair()
        #print i1, i2
        j1 = self.replicas[i1]
        j2 = self.replicas[i2]
        e1 = self.energies[j1]
        e2 = self.energies[j2]
        beta1 = self.beta[i1]
        beta2 = self.beta[i2]

        w = min( 1. , np.exp( (e1-e2) * (beta1-beta2) ) )
        rand = np.random.rand()
        if w > rand:
            #accept step
            #print "accept step"
            k = self.replicas[i1]
            self.replicas[i1] = self.replicas[i2]
            self.replicas[i2] = k
            self.naccept += 1
            self.acccount[[i1,i2]] += 1
        else:
            #reject step
            #print "reject step"
            pass