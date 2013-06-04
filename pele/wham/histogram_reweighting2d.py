from numpy import exp, log
import numpy as np #to access np.exp() not built int exp
#from math import *
#import scipy.optimize
#import timeseries # for timeseries analysis 
import wham_utils
from wham_utils import logSum, logSum1
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
#mbar = pickle.load(open("mbar.pickle","rb"))






class wham2d(object):
    """ class to combine 2d histograms of energy E and order parameter q at
    multiple temperatures into one best estimate for the histogram
  
    input will be:
    filenames : list of filenames where the data can be found
    Tlist # Tlist[k] is the temperature of the simulation in filenames[k] 
  
    binenergy = zeros(nebins, float64) #lower edge of bin energy
    binq =      zeros(nqbins, float64) #lower edge of q bins
    visits2d =  zeros([nrep,nebins,nqbins], integer) #2d histogram of data
    """
    #=============================================================================================
    # Constructor.
    #=============================================================================================
    def __init__(self, Tlist, binenergy, binq, visits2d):
    
        #define some parameters
        self.k_B=1.
        self.LOGMIN = -1e200
    
        self.nrep = len(Tlist)
        self.nebins = len(binenergy)
        self.nqbins = len(binq)
        self.Tlist = np.array(Tlist)
        self.binenergy = np.array(binenergy)
        self.binq = np.array(binq)
        self.visits2d = np.array(visits2d, dtype = np.integer)
        
    def minimize(self):
        #shape(visits2d) is now (nqbins, nebins, nreps)
        #we need it to be (nreps, nqbins*nebins)
        #first reorder indices
        nreps = self.nrep
        nebins = self.nebins
        nqbins = self.nqbins
        nbins = self.nebins * self.nqbins
        #visits = np.zeros([nreps, nebins, nqbins], np.integer)
        reduced_energy = np.zeros([nreps, nebins, nqbins])
#        for k in range(self.nrep):
#            for j in range(self.nqbins):
#                for i in range(self.nebins):
#                    #visits[k,i,j] = self.visits2d[i,j,k]
#                    reduced_energy[k,i,j] = self.binenergy[i] / (self.Tlist[k]*self.k_B)
        for j in range(self.nqbins):
            reduced_energy[:,:,j] = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis]*self.k_B)
                    
        visits = self.visits2d
        visits = np.reshape( visits, [nreps, nbins ]) 
        reduced_energy = np.reshape( reduced_energy, [nreps, nbins])           
        self.logP = np.where( visits != 0, np.log( visits.astype(float) ), 0 )

        
        from wham_potential import WhamPotential
        whampot = WhamPotential( self.logP, reduced_energy )
        
        nvar = nbins + nreps
        X = np.random.rand(nvar)
        print "initial energy", whampot.getEnergy(X)
        try: 
            from pele.optimize import mylbfgs as quench
            ret = quench(X, whampot, iprint=10, maxstep = 1e4)
        except ImportError:
            from pele.optimize import lbfgs_scipy as quench
            ret = quench(X, whampot)            

        print "quenched energy", ret.energy
        
        global_min = False
        if global_min:
            from pele.basinhopping import BasinHopping
            from pele.takestep.displace import RandomDisplacement
            takestep = RandomDisplacement(stepsize=10.)
            takestep.useAdaptiveStep()
            takestep.adaptive_class.f = 2.
            bh = BasinHopping(X, whampot, takestep)
            bh.run(1000)
        
        
        #self.logn_Eq = zeros([nebins,nqbins], float64)
        X = ret.coords
        self.logn_Eq = X[nreps:]
        self.w_i_final = X[:nreps]
        
        self.logn_Eq = np.reshape(self.logn_Eq, [nebins, nqbins])
        self.logn_Eq = np.where( self.visits2d.sum(0) == 0, self.LOGMIN, self.logn_Eq )
        
        #renormalize logn_Eq
        #self.allzero2dind = np.where(self.visits2d.sum(2) == 0)
        #self.notzero2dind = np.where(self.visits2d.sum(2) != 0)
        #print self.allzero2dind
        #self.logn_Eq -= np.min(self.logn_Eq[ self.notzero2dind[1], self.notzero2dind[0] ])

        #self.logn_Eq[self.allzero2dind] = 0



  

  
    def calc_Fq(self, TRANGE = []):
        self.allzero2dind = np.where(self.visits2d.sum(0) == 0)

  
        #put some variables in this namespace
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        if len(TRANGE) == 0:
            NTEMP = 5 # number of temperatures to calculate expectation values
            TMAX = self.Tlist[-1]
            TMIN = self.Tlist[0]
            TINT=(TMAX-TMIN)/(NTEMP-1)
            TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
    
        #find the occupied bin with the minimum energy
        EREF=0
        for i in range(nebins):
            if visits2d[:,i,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        self.nodataq = np.where((visits2d.sum(0).sum(0)) == 0)
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        self.F_q = np.zeros([nqbins,len(TRANGE)])
        F_q = self.F_q
        logP_Eq = np.zeros([nebins,nqbins])
        logP_q = np.zeros(nqbins)
        for n in range(len(TRANGE)):
            T=TRANGE[n]
            for i in range(nebins):
                logP_Eq[i,:] = logn_Eq[i,:]-(binenergy[i] - EREF)/(self.k_B*T)
      
            logP_Eq[self.allzero2dind[0], self.allzero2dind[1]] = self.LOGMIN
            expoffset = np.nanmax(logP_Eq)
            #print "T expoffset ", T, expoffset
            logP_Eq -= expoffset
            #P_q = np.exp(logP_Eq).sum(0)
            # sum over the energy
            for j in range(nqbins):
                logP_q[j] = logSum( logP_Eq[:,j] )
            logP_q[self.nodataq] = np.NaN
            F_q[:,n] = -self.k_B*T*logP_q[:]
            fmin = np.nanmin(F_q[:,n])
            F_q[:,n] -= fmin
    
        return TRANGE,F_q
  
    def calc_qavg(self, TRANGE = []):
        """calculate the average q as a function of temperature"""
        #put some variables in this namespace
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        binq=self.binq
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        if len(TRANGE) == 0:
            NTEMP = 100 # number of temperatures to calculate expectation values
            TMAX = self.Tlist[-1]
            TMIN = self.Tlist[0]
            TINT=(TMAX-TMIN)/(NTEMP-1)
            TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
    
        #find the ocupied bin with the minimum energy
        EREF=0
        for i in range(nebins):
            if visits2d[:,i,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        #don't need to recalculate it
        #self.nodataq = where((visits2d.sum(2).sum(0)) == 0)
    
        #calculate the mean q at each temperature
        self.qavg = np.zeros(len(TRANGE))
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        logP_Eq = np.zeros([nebins,nqbins])
        logP_q = np.zeros(nqbins)
        for n in range(len(TRANGE)):
            T=TRANGE[n]
            for i in range(nebins):
                logP_Eq[i,:] = logn_Eq[i,:]-(binenergy[i] - EREF)/(self.k_B*T)
      
            logP_Eq[self.allzero2dind[0], self.allzero2dind[1]] = self.LOGMIN
            expoffset = logP_Eq.max()
            #print "T expoffset ", T, expoffset
            logP_Eq -= expoffset
            #P_q = np.exp(logP_Eq).sum(0)
            # sum over the energy
            for j in range(nqbins):
                logP_q[j] = wham_utils.logSum( logP_Eq[:,j] )
            logP_q[self.nodataq] = np.NaN
      
            #find mean q
            qmin = min(binq)
            qmin -= 0.1
            lqavg = -1.0e30
            lnorm = -1.0e30
            for i in range(0,nqbins): 
                if not np.isnan(logP_q[i]):
                    lnorm = wham_utils.logSum1( lnorm, logP_q[i] ) 
                    lqavg = wham_utils.logSum1( lqavg, logP_q[i] + log(binq[i] - qmin) )
            self.qavg[n] = exp(lqavg - lnorm) + qmin
            #print lqavg
    
        return TRANGE,self.qavg
  
  
    def calc_Cv(self, NDOF):
        nebins = self.nebins
        visits1d = self.visits2d.sum(2)
        logn_E = np.zeros(nebins)
        for i in range(nebins):
            logn_E[i] = wham_utils.logSum(self.logn_Eq[i,:])

        return wham_utils.calc_Cv(logn_E, visits1d, self.binenergy, \
                                  NDOF, self.Tlist, self.k_B)

  

