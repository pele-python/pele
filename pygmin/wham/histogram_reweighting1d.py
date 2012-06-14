import numpy as np #to access np.exp() not built int exp
#import timeseries # for timeseries analysis 
#import commands
#import pdb;
#import pickle
from wham_potential import WhamPotential
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import *
import wham_utils






class wham1d:
    """ class to combine 1d histograms of energy E at
    multiple temperatures into one best estimate for the histogram

    input will be:
    filenames : list of filenames where the data can be found
    Tlist # Tlist[k] is the temperature of the simulation in filenames[k] 

    binenergy = zeros(nebins, float64) #lower edge of bin energy
    visits1d =  zeros([nebins,nrep], integer) #1d histograms of data
    """
    #=============================================================================================
    # Constructor.
    #=============================================================================================
    def __init__(self, Tlist, binenergy, visits1d):

        #define some parameters
        self.k_B=1.

        self.nrep = len(Tlist)
        self.nebins = len(binenergy)
        self.Tlist = np.array(Tlist, dtype = np.float64)
        self.binenergy = np.array(binenergy, dtype = np.float64)
        self.visits1d = np.array(visits1d, dtype = np.integer)

    def minimize(self):
        nreps = self.nrep
        nbins = self.nebins
        visitsT = (self.visits1d)
        #print "min vis", np.min(visitsT)
        self.logP = np.where( visitsT != 0, np.log( visitsT ), 0 )
        #print "minlogp", np.min(self.logP)
        self.reduced_energy = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis] * self.k_B)
        
        self.whampot = WhamPotential(self.logP, self.reduced_energy)
        
        X = np.random.rand( nreps + nbins )
        E = self.whampot.getEnergy(X)
        print "energy", E 
        
        print "quenching"
        from pygmin.optimize.quench import quench
        ret = quench(X, self.whampot.getEnergyGradient)
        print "quench energy", ret[1]
        X = ret[0]
        
        self.logn_E = X[nreps:]
        self.w_i_final = X[:nreps]
        
    def globalMinimization(self):
        """
        in experimentation i've never been able to find more than
        one minimum
        """
        nreps = self.nrep
        nbins = self.nebins
        visitsT = (self.visits1d)
        #print "min vis", np.min(visitsT)
        self.logP = np.where( visitsT != 0, np.log( visitsT ), 0 )
        #print "minlogp", np.min(self.logP)
        self.reduced_energy = self.binenergy[np.newaxis,:] / (self.Tlist[:,np.newaxis] * self.k_B)
        
        self.whampot = WhamPotential(self.logP, self.reduced_energy)
        
        X = np.random.rand( nreps + nbins )
        E = self.whampot.getEnergy(X)
        print "energy", E 
        
        print "quenching"
        from pygmin.optimize.quench import quench
        ret = quench(X, self.whampot.getEnergyGradient)
        print "quench energy", ret[1]
        
        from pygmin.basinhopping import BasinHopping
        from pygmin.takestep.displace import RandomDisplacement
        takestep = RandomDisplacement(stepsize=10)
        takestep.useAdaptiveStep()
        takestep.adaptive_class.f = 1.5 #i have no idea what a good stepsize should be
        bh = BasinHopping(X, self.whampot, takestep )
        
        import matplotlib.pyplot as plt
        for i in range(10):
            bh.run(2000)
            self.logn_E = bh.coords[nreps:]
            cvdata = self.calc_Cv(400)
            plt.plot(cvdata[:,0], cvdata[:,5], '-')
        plt.show()
            
            
        
        X = bh.coords
        self.logn_E = X[nreps:]
        self.w_i_final = X[:nreps]




    def calc_Cv(self, NDOF, TRANGE=None, NTEMP=100):
        return wham_utils.calc_Cv(self.logn_E, self.visits1d, self.binenergy, \
                NDOF, self.Tlist, self.k_B, TRANGE, NTEMP)



