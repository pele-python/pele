import numpy as np


class manageStepSize:
    """a class to manage the adaptive step size"""
    def __init__(self, stepsize, accrat, nstepsaccrat=100, f = 1.05):
        self.stepsize = stepsize
        self.f = f
        self.accrat = accrat #target accept ratio
        self.nstepsaccrat = nstepsaccrat #interval for adjusting septsize

        self.nsteps = 0
        self.nstepstot = 0
        self.naccepted = 0
        self.naccepted_tot = 0
        self.nadjust = 0
        self.changehist = []

    def getStepSize(self):
        return self.stepsize

    def adjustFactor(self, factor):
        """This is the function which ensures f is a reasonable number.  Normal
        GMIN doesn't adjust f.  So this function currently does nothing
        """
        return 
        changef = 1.5
        if len(self.changehist) == 0: return 
        if (factor < 1.0) != (self.changehist[-1] < 1.0):
            self.f -= (self.f - 1)/1.5


    def adjustStep(self):
        """adjust stepsize"""
        self.nadjust += 1
        rat = float(self.naccepted)/self.nsteps
        if rat < self.accrat:
            #reduce step size
            self.adjustFactor( 1./self.f )
            self.stepsize /= self.f
            self.changehist.append( 1./self.f )
        else:
            #increase step size
            self.adjustFactor( self.f )
            self.stepsize *= self.f
            self.changehist.append( self.f )
        self.nsteps = 0
        self.nsteps = 0
        self.naccepted = 0
        print "accrat was ", rat, "new stepsize is ", self.stepsize, "f is", self.f

    def insertStep(self, accepted ):
        """tell us whether a step was accepted or rejected"""
        if accepted: 
            self.naccepted += 1
            self.naccepted_tot += 1
        self.nsteps += 1
        self.nstepstot += 1
        if self.nsteps == self.nstepsaccrat:
            self.adjustStep()

    def insertStepWrapper(self, E, coords, accepted ):
        return self.insertStep(accepted)

    def accratTot(self):
        return float(self.naccepted_tot) / self.nstepstot
