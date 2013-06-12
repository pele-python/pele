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
    visits1d =  zeros([nrep,nebins], integer) #1d histograms of data
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
        #print "energy", E 
        
        #print "quenching"
        try: 
            from pele.optimize import mylbfgs as quench
            ret = quench(X, self.whampot, iprint=-1, maxstep=1e4)
        except ImportError:
            from pele.optimize import lbfgs_scipy as quench
            ret = quench(X, self.whampot)            
        #print "quench energy", ret.energy
        X = ret.coords
        
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
        from pele.optimize import lbfgs_scipy as quench
        ret = quench(X, self.whampot)
        print "quench energy", ret.energy
        
        from pele.basinhopping import BasinHopping
        from pele.takestep.displace import RandomDisplacement
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

#    def calc_Cv_no_wham(self):
#        """ """

    def calc_Cv_new(self, NDOF, TRANGE=None, NTEMP=100):
        from pele.thermodynamics import dos_to_cv
        dT = (self.Tlist[-1] - self.Tlist[0]) / NTEMP
        Tlist = np.arange(self.Tlist[0], self.Tlist[-1], dT)
#        print self.logn_E
        lZ, U, U2, Cv = dos_to_cv(self.binenergy, self.logn_E, Tlist, K=NDOF)
        cvdata = np.zeros([len(Tlist), 6])
        cvdata[:,0] = Tlist
        cvdata[:,1] = lZ
        cvdata[:,2] = U # average potential energy
        cvdata[:,3] = U2
        cvdata[:,5] = Cv
        
        eavg = U + float(NDOF) / 2 * Tlist # average energy including the kinetic degrees of freedom
        cvdata[:,4] = eavg
        return cvdata
        

    def calc_Cv(self, NDOF, TRANGE=None, NTEMP=100, use_log_sum=None):
#        return self.calc_Cv_new(NDOF, TRANGE, NTEMP)
        return wham_utils.calc_Cv(self.logn_E, self.visits1d, self.binenergy,
                NDOF, self.Tlist, self.k_B, TRANGE, NTEMP, use_log_sum=use_log_sum)



