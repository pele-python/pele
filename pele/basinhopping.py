# -*- coding: iso-8859-1 -*-
import sys
from pele.mc import MonteCarlo
from pele.optimize import mylbfgs

class BasinHopping(MonteCarlo):
    """
    A class to run the basin hopping algorithm

    Parameters
    ----------
    All required and optional parameters from base class MonteCarlo :
    quench : callable, optional
        Use this quencher as default
    insert_rejected : bool
        insert the rejected structure into the storage class
    
    Notes
    -----
    Basin-hopping is a stochastic algorithm which attempts to find the global
    minimum of a smooth scalar function of one or more variables

    The algorithm is iterative with each cycle composed of the following
    features

    1) random perturbation of the coordinates

    2) local minimization

    3) accept or reject the new coordinates based on the minimized function
       value

    The acceptance test used here is the Metropolis criterion of standard Monte
    Carlo algorithms, although there are many other possibilities [3]_.

    This global minimization method has been shown to be extremely efficient
    for a wide variety of problems in physics and chemistry.  It is
    particularly useful when the function has many minima separated by large
    barriers. See the Cambridge Cluster Database
    http://www-wales.ch.cam.ac.uk/CCD.html for databases of molecular systems
    that have been optimized primarily using basin-hopping.  This database
    includes minimization problems exceeding 300 degrees of freedom.
    
    See Also
    --------
    pele.mc.MonteCarlo : base class
    pele.potentials, pele.takestep, pele.storage, pele.accept_tests


    """

    def __init__(self, coords, potential, takeStep, storage=None, event_after_step=None, acceptTest=None,
                 temperature=1.0, quench=None, confCheck=None, outstream=sys.stdout, insert_rejected=False):
        #########################################################################
        # initialize MonteCarlo base class
        #########################################################################
        MonteCarlo.__init__(self, coords, potential, takeStep, storage=storage,
                            event_after_step=event_after_step,
                            acceptTest=acceptTest, temperature=temperature,
                            confCheck=confCheck, outstream=outstream,
                            store_initial=False)

        if quench is None:
            quench = lambda coords : mylbfgs(coords, self.potential)
        self.quench = quench
                
        #########################################################################
        # do initial quench
        #########################################################################
        self.markovE_old = self.markovE
        res = self.quench(self.coords)
        self.result.nfev += res.nfev

        self.coords = res.coords
        self.markovE = res.energy
        self.rms = res.rms
        self.funcalls = res.nfev

        self.insert_rejected = insert_rejected
        
        if self.storage:
            self.storage(self.markovE, self.coords)
        
        # print the initial quench
        self.acceptstep = True
        self.trial_energy = self.markovE
        self.printStep()
        
        self.result.energy = self.markovE
        self.result.coords = self.coords.copy()


    def _mcStep(self):
        """
        take one monte carlo basin hopping step

        overload the MonteCarlo base class step
        """
        self.coords_after_step = self.coords.copy() # make  a working copy
        #########################################################################
        # take step
        #########################################################################
        self.takeStep.takeStep(self.coords_after_step, driver=self)

        #########################################################################
        # quench
        #########################################################################
        res = self.quench(self.coords_after_step)
        self.result.nfev += res.nfev
#        if isinstance(res, tuple): # for compatability with old and new quenchers
#            res = res[4]
        self.trial_coords = res.coords
        self.trial_energy = res.energy
        self.rms = res.rms
        self.funcalls = res.nfev

        #########################################################################
        # check if step is a valid configuration, otherwise reject
        #########################################################################
        self.acceptstep = True
        self.config_ok = True
        for check in self.confCheck:
            if not check(self.trial_energy, self.trial_coords, driver=self):
                self.acceptstep=False
                self.config_ok = False
        
        #########################################################################
        # check whether step is accepted with user defined tests.  If any returns
        # false then reject step.
        #########################################################################
        if self.acceptstep:
            self.acceptstep = self.acceptTest(self.markovE, self.trial_energy,
                                              self.coords, self.trial_coords)

        #########################################################################
        # return new coords and energy and whether or not they were accepted
        #########################################################################
        return self.acceptstep, self.trial_coords, self.trial_energy


    def printStep(self):
        if self.stepnum % self.printfrq == 0:
            if self.outstream != None:
                self.outstream.write("Qu   " + str(self.stepnum) + " E= " +
                                     str(self.trial_energy) + " quench_steps= " +
                                     str(self.funcalls) + " RMS= " + str(self.rms) + 
                                     " Markov E= " + str(self.markovE_old) + 
                                     " accepted= " + str(self.acceptstep) + "\n" )
    
    def __getstate__(self):
        ddict = self.__dict__.copy()
        del ddict["outstream"]
        del ddict["potential"]
        return ddict
    
    def __setstate__(self, dct):
        self.__dict__.update(dct)
        self.outstream = sys.stdout


if __name__ == "__main__":
    from pele.systems import LJCluster
    natoms = 13
    sys = LJCluster(natoms)
    bh = sys.get_basinhopping()
    bh.run(100)

