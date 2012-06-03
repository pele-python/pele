from numpy import *
import numpy as np #to access np.exp() not built int exp
from math import *
import scipy.optimize
#import timeseries # for timeseries analysis 
import commands
import pdb;
import pickle
import load_data
import wham_utils
from wham_utils import logSum, logSum1
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
#mbar = pickle.load(open("mbar.pickle","rb"))

class loadData2dExp:
    def __init__(self, filenames, ecolumn, qcolumn, nqbins, qmin = 0, qmax = 1.0, nebins = 200, NEGLECT = 0.01, fskip=0., qcombine=[], dEmin=0.1):

        nrep = len(filenames)

        #if the columns have a header, print the name of the column requested
        with open(filenames[0], "r") as fin:
            line = fin.readline()
            sline = line.split()
            if sline[0][0] == '#':
                print "qcolumn header ", sline[qcolumn]


        #qcombine = [qcolumn1, qcolumn2, frac]
        #make a linear combination of order parameters from columns 1 and 2
        # q = (q1 + q2weight*q2)/(1+q2weight)
        qcomb = False
        if len(qcombine) == 3:
            qcomb = True
            qcolumn1 = int(qcombine[0])
            qcolumn2 = int(qcombine[1])
            q2weight = qcombine[2]   
            print "qcombine: ", qcolumn1, qcolumn2, q2weight

        #loop through files and get maximum and minimum energies
        print "determining emax and emin from input data"
        emin=1e10
        emax=-1e10
        nlines=[]
        Ebar=[]
        Esig=[]
        for f in filenames:
            data=genfromtxt(f)
            e=data[:,ecolumn].min()
            if e < emin: emin = e
            e=data[:,ecolumn].max()
            if e > emax: emax = e
            nlines.append(len(data[:,0]))
            ebar = mean( data[:,ecolumn] )
            esig = std( data[:,ecolumn] )
            Ebar.append(ebar)
            Esig.append(esig)

        emax += 1e-4 #so the highest energy returns the last bin, not one past the last bin
        qmax += 1e-4 #so the highest q returns the last bin, not one past the last bin

        dE = (emax-emin)/(nebins)
        dq = (qmax-qmin)/(nqbins)
        print emin, " < E < ", emax, " dE = ", dE

        self.visits2d = zeros([nebins,nqbins,nrep], int32) 
        self.binenergy = array([ emin + dE*i for i in range(nebins)])
        self.binq = array([ qmin + dq*i for i in range(nqbins)])
        visits2d = self.visits2d
        binenergy = self.binenergy
        binq = self.binq

        #for histogram, numpy wants bin endges, including rightmost edge
        binenergynp = array([ emin + dE*i for i in range(nebins+1)]) 
        binqnp = array([ qmin + dq*i for i in range(nqbins+1)])


        #use exponentially increasing bin width
        exponential_bins = True
        if exponential_bins:
            load_data.exponentialBinEnergy( binenergynp, emin, emax, nebins, Ebar[0], Esig[0], Ebar[-1], Esig[-1] )
            binenergy[0:nebins] = ( binenergynp[0:-1])

        #load the energies into visits2d
        print "reading data from files"
        for k in range(nrep):
            with open(filenames[k],"r") as f:
                data=genfromtxt(f)
                lfirst = nlines[k]*fskip  
                e = data[:,ecolumn]
                if qcomb:
                    q1 = float(line.split()[qcolumn1])
                    q2 = float(line.split()[qcolumn2])
                    q = (float(q1) + q2weight*q2 )/(1.+q2weight)
                    q = (data[:,qcolumn1] + data[:,qcolumn2]*q2weight) / (1+q2weight)
                else:
                    q = data[:,qcolumn]
                visits2d[:,:,k], xbins, ybins = np.histogram2d(e, q, [binenergynp, binqnp])

        #set visits to zero if the count is below a threshold
        #the threshold is different for each replica, and is determined by the ratio to the maximum count
        for k in range(nrep):
            visitsmax=visits2d[:,:,k].max()
            minvis = max(1, int(visitsmax * NEGLECT))
            ind=where(visits2d[:,:,k] <= minvis)
            print k, "total visits ", visits2d[:,:,k].sum(), "max visits ", visitsmax, "max visits to be discarded ", visits2d[ind[0],ind[1],k].max()
            visits2d[ind[0],ind[1],k] = 0





class wham2d:
    """ class to combine 2d histograms of energy E and order parameter q at
    multiple temperatures into one best estimate for the histogram
  
    input will be:
    filenames : list of filenames where the data can be found
    Tlist # Tlist[k] is the temperature of the simulation in filenames[k] 
  
    binenergy = zeros(nebins, float64) #lower edge of bin energy
    binq =      zeros(nqbins, float64) #lower edge of q bins
    visits2d =  zeros([nebins,nqbins,nrep], integer) #2d histogram of data
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
        self.Tlist = array(Tlist, dtype = float64)
        self.binenergy = array(binenergy, dtype = float64)
        self.binq = array(binq, dtype = float64)
        self.visits2d = array(visits2d, dtype = int32)
        
    def minimize(self):
        #shape(visits2d) is now (nqbins, nebins, nreps)
        #we need it to be (nreps, nqbins*nebins)
        #first reorder indices
        nreps = self.nrep
        nebins = self.nebins
        nqbins = self.nqbins
        nbins = self.nebins * self.nqbins
        visits = np.zeros([nreps, nebins, nqbins], np.integer)
        reduced_energy = np.zeros([nreps, nebins, nqbins])
        for k in range(self.nrep):
            for j in range(self.nqbins):
                for i in range(self.nebins):
                    visits[k,i,j] = self.visits2d[i,j,k]
                    reduced_energy[k,i,j] = self.binenergy[i] / (self.Tlist[k]*self.k_B)
                    
        
        visits = np.reshape( visits, [nreps, nbins ]) 
        reduced_energy = np.reshape( reduced_energy, [nreps, nbins])           
        self.logP = np.where( visits != 0, np.log( visits.astype(float) ), 0 )
        #assert((visits >= 0).all())
        #print visits
        #print np.log(np.max(visits))
        #print self.logP
        #print np.max(visits), np.max(self.visits2d)
        
        from optimize.quench import quench
        from wham_potential import WhamPotential
        whampot = WhamPotential( self.logP, reduced_energy )
        
        nvar = nbins + nreps
        X = np.random.rand(nvar)
        print "initial energy", whampot.getEnergy(X)
        ret = quench(X, whampot.getEnergyGradient)
        print "quenched energy", ret[1]
        
        global_min = False
        if global_min:
            from basinhopping import BasinHopping
            from take_step.random_displacement import takeStep
            takestep = takeStep(stepsize=10.)
            takestep.useAdaptiveStep()
            takestep.adaptive_class.f = 2.
            bh = BasinHopping(X, whampot, takestep)
            bh.run(1000)
        
        
        #self.logn_Eq = zeros([nebins,nqbins], float64)
        X = ret[0]
        self.logn_Eq = X[nreps:]
        self.w_i_final = X[:nreps]
        
        self.logn_Eq = np.reshape(self.logn_Eq, [nebins, nqbins])
        
        #renormalize logn_Eq
        #self.logn_Eq -= np.min(self.logn_Eq)
        #self.allzero2dind = where(self.visits2d.sum(2) == 0)

        #self.logn_Eq[self.allzero2dind] = 0



  

  
    def calc_Fq(self, TRANGE = []):
        self.allzero2dind = where(self.visits2d.sum(2) == 0)

  
        #put some variables in this namespace
        nrep=self.nrep
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        binq=self.binq
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        if len(TRANGE) == 0:
            NTEMP = 5 # number of temperatures to calculate expectation values
            TMAX = self.Tlist[-1]
            TMIN = self.Tlist[0]
            TINT=(TMAX-TMIN)/(NTEMP-1)
            TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
    
        #find the ocupied bin with the minimum energy
        EREF=0
        for i in range(nebins):
            if visits2d[i,:,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        self.nodataq = where((visits2d.sum(2).sum(0)) == 0)
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        self.F_q = zeros([nqbins,len(TRANGE)], float64)
        F_q = self.F_q
        logP_Eq = zeros([nebins,nqbins], float64)
        logP_q = zeros(nqbins, float64)
        for n in range(len(TRANGE)):
            T=TRANGE[n]
            for i in range(nebins):
                logP_Eq[i,:] = logn_Eq[i,:]-(binenergy[i] - EREF)/(self.k_B*T)
      
            logP_Eq[self.allzero2dind[0], self.allzero2dind[1]] = self.LOGMIN
            expoffset = nanmax(logP_Eq)
            print "T expoffset ", T, expoffset
            logP_Eq -= expoffset
            #P_q = np.exp(logP_Eq).sum(0)
            # sum over the energy
            for j in range(nqbins):
                logP_q[j] = logSum( logP_Eq[:,j] )
            logP_q[self.nodataq] = NaN
            F_q[:,n] = -self.k_B*T*logP_q[:]
            fmin = nanmin(F_q[:,n])
            F_q[:,n] -= fmin
    
        return TRANGE,F_q
  
    def calc_qavg(self, TRANGE = []):
        """calculate the average q as a function of temperature"""
        #put some variables in this namespace
        nrep=self.nrep
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
            if visits2d[i,:,:].sum() > 0:
                EREF = binenergy[i]
                break
    
        #don't need to recalculate it
        #self.nodataq = where((visits2d.sum(2).sum(0)) == 0)
    
        #calculate the mean q at each temperature
        self.qavg = zeros(len(TRANGE), float64)
    
        #now calculate P(q,T)
        # P(q,T) = sum_E n(E,q)*exp(-E/T)  
        #TRANGE=range(1,9)
        logP_Eq = zeros([nebins,nqbins], float64)
        logP_q = zeros(nqbins, float64)
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
            logP_q[self.nodataq] = NaN
      
            #find mean q
            qmin = min(binq)
            qmin -= 0.1
            lqavg = -1.0e30
            lnorm = -1.0e30
            for i in range(0,nqbins): 
                if not isnan(logP_q[i]):
                    lnorm = wham_utils.logSum1( lnorm, logP_q[i] ) 
                    lqavg = wham_utils.logSum1( lqavg, logP_q[i] + log(binq[i] - qmin) )
            self.qavg[n] = exp(lqavg - lnorm) + qmin
            #print lqavg
    
        return TRANGE,self.qavg
  
  
    def calc_Cv(self, NDOF, fout):
  
        #put some variables in this namespace
        nrep=self.nrep
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        binq=self.binq
        visits2d=self.visits2d
        logn_Eq=self.logn_Eq
    
        logn_E = zeros(nebins, float64)
        for i in range(nebins):
            logn_E[i] = wham_utils.logSum(logn_Eq[i,:])
        self.nodatae = where((visits2d.sum(2).sum(1)) == 0)
        self.allzeroe = (visits2d.sum(2).sum(1)) == 0
        #fout=open("weights.A2d","w")
        #savetxt(fout,column_stack((binenergy,logn_E)))
        #fout.close()
    
        #find the ocupied bin with the minimum energy
        for i in range(nebins):
            if not self.allzeroe[i] :
                EREF = binenergy[i]
                break
    
        #now calculate partition functions, energy expectation values, and Cv
        #fout = open("Cv.out.Apy2d", "w")
        #NDOF = 12*3 # = number of degrees of freedom
        NTEMP = 100 # number of temperatures to calculate expectation values
        TMAX = self.Tlist[-1]
        TMIN = self.Tlist[0]
        TINT=(TMAX-TMIN)/(NTEMP-1)
        TRANGE = [ TMIN + i*TINT for i in range(NTEMP) ]
        for T in TRANGE:
            kBT = self.k_B*T
            Z0=0.0
            Z1=0.0
            Z2=0.0
            #find expoffset so the exponentials don't blow up
            expoffset=-1e10
            for i in range(nebins):
                if self.allzeroe[i]: continue
                EDIFF = (binenergy[i]-EREF)
                dummy = ( logn_E[i]   -(EDIFF)/(kBT))
                if dummy > expoffset: expoffset = dummy
            #do calculation
            for i in range(nebins):
                if self.allzeroe[i]: continue
                EDIFF = (binenergy[i]-EREF)
                dummy = exp( logn_E[i]   -(EDIFF)/(kBT) - expoffset)
                Z0 += dummy
                Z1 += dummy * EDIFF
                Z2 += dummy * EDIFF * EDIFF
      
            if i == nebins-1:
                dE = binenergy[i]-binenergy[i-1]
            else:
                dE = binenergy[i+1]-binenergy[i]
            if dE/kBT < 1.0E-7:
                ONEMEXP=-dE/kBT
            else:
                ONEMEXP= 1.0-exp(dE/kBT)
      
            Eavg = NDOF*kBT/2.0 + 1.0*(kBT + dE/ONEMEXP) + Z1/Z0 + EREF
            Cv = NDOF/2. + 1.0*(1.0 - dE**2 * exp(dE/kBT)/(ONEMEXP**2*kBT**2)) - (Z1/(Z0*kBT))**2 + Z2/(Z0*kBT**2)
            array([kBT, Z0*exp(expoffset), Z1*exp(expoffset), Z2*exp(expoffset), Eavg, Cv, log(Z0)+expoffset, log(Z1)+expoffset, log(Z2)+expoffset]).tofile(fout," ")
            fout.write("\n")
  

