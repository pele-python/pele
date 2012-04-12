import numpy as np

class Metropolis:
    def __init__(self, temperature, random=np.random.rand):
        self.random = random
        self.temperature = temperature

    def acceptReject(self, Eold, Enew, newcoords=[]):
        acceptstep = True
        wcomp = (Enew - Eold)/self.temperature
        w=min(1.0,np.exp(-wcomp))
        rand = self.random()
        if (rand > w): acceptstep = False

        #print "mc step: Eo", Eold, "En", Enew, "accepted", acceptstep

        return acceptstep
        
