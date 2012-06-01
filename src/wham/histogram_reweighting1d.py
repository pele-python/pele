import numpy as np #to access np.exp() not built int exp
#import timeseries # for timeseries analysis 
#import commands
#import pdb;
#import pickle
from wham_potential import WhamPotential
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *






class wham1d:
    """ class to combine 1d histograms of energy E at
    multiple temperatures into one best estimate for the histogram

    input will be:
    filenames : list of filenames where the data can be found
    Tlist # Tlist[k] is the temperature of the simulation in filenames[k] 

    binenergy = zeros(nebins, float64) #lower edge of bin energy
    visits1d =  zeros([nebins,nqbins,nrep], integer) #1d histogram of data
    """
    #=============================================================================================
    # Constructor.
    #=============================================================================================
    def __init__(self, Tlist, binenergy, visits1d):

        #define some parameters
        self.k_B=1.
        self.LOGMIN = -1e200

        self.nrep = len(Tlist)
        self.nebins = len(binenergy)
        self.nqbins = 1
        self.Tlist = np.array(Tlist, dtype = np.float64)
        self.binenergy = np.array(binenergy, dtype = np.float64)
        self.visits1d = np.array(visits1d, dtype = np.int32)

    def minimizeNew(self):
        nreps = self.nrep
        nbins = self.nebins
        visitsT = np.transpose(self.visits1d)
        #print "min vis", np.min(visitsT)
        self.logP = np.where( visitsT != 0, np.log( visitsT ), 0 )
        #print "minlogp", np.min(self.logP)
        self.reduced_energy = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis] * self.k_B)
        
        self.whampot = WhamPotential(self.logP, self.reduced_energy)
        
        X = np.random.rand( nreps + nbins )
        E = self.whampot.getEnergy(X)
        print "energy", E 
        
        print "quenching"
        from optimize.quench import quench
        ret = quench(X, self.whampot.getEnergyGradient)
        print "quench energy", ret[1]
        X = ret[0]
        self.logn_E = X[nreps:]
        self.w_i_final = X[:nreps]
        
        


    def calc_Cv(self, NDOF, fout):
        #put some variables in this namespace
        nrep=self.nrep
        nebins=self.nebins
        nqbins=self.nqbins
        binenergy=self.binenergy
        visits1d=self.visits1d
        logn_E=self.logn_E

        self.nodatae = np.where((visits1d.sum(1)) == 0)
        self.allzeroe = (visits1d.sum(1)) == 0
        #fout=open("weights.A1d","w")
        #savetxt(fout,column_stack((binenergy,logn_E)))
        #fout.close()

        #find the ocupied bin with the minimum energy
        for i in range(nebins):
            if not self.allzeroe[i] :
                EREF = binenergy[i]
                break

        #now calculate partition functions, energy expectation values, and Cv
        #fout = open("Cv.out.Apy1d", "w")
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
                dummy = np.exp( logn_E[i]   -(EDIFF)/(kBT) - expoffset)
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
                ONEMEXP= 1.0-np.exp(dE/kBT)

            #Eavg = NDOF*kBT/2.0 + 1.0*(kBT + dE/ONEMEXP) + Z1/Z0 + EDIFF
            Eavg = NDOF*kBT/2.0 + 1.0*(kBT + dE/ONEMEXP) + Z1/Z0 + EREF
            Cv = NDOF/2. + 1.0*(1.0 - dE**2 * np.exp(dE/kBT)/(ONEMEXP**2*kBT**2)) - (Z1/(Z0*kBT))**2 + Z2/(Z0*kBT**2)
            np.array([kBT, Z0*np.exp(expoffset), Z1*np.exp(expoffset), Z2*np.exp(expoffset), Eavg, Cv, np.log(Z0)+expoffset, np.log(Z1)+expoffset, np.log(Z2)+expoffset]).tofile(fout," ")
            fout.write("\n")


