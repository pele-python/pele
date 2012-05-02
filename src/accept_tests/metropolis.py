import numpy as np

class Metropolis:
    def __init__(self, temperature, random=np.random.rand):
        self.random = random
        self.temperature = temperature

    def acceptRejectE(self, Eold, Enew):
        acceptstep = True
        wcomp = (Enew - Eold)/self.temperature
        w=min(1.0,np.exp(-wcomp))
        rand = self.random()
        if (rand > w): acceptstep = False

        #print "mc step: Eo", Eold, "En", Enew, "accepted", acceptstep

        return acceptstep

    def acceptReject(self, Eold, Enew, qcoords=[], coords=[]):
        return self.acceptRejectE(Eold, Enew)

class MetropolisNonQuench:
    """
    perform metropolis criterion on non quenched energy
    """
    def __init__(self, temperature, potential, random=np.random.rand):
        self.potential = potential
        self.metropolis = Metropolis(temperature, random)


    def acceptReject(self, Equench_old=0., Equench_new=0., qcoords=[], coords=[]):
        self.Enew = self.potential.getEnergy(coords)
        try:
            self.Eold
        except:
            self.Eold = self.Enew
        #print "Eold Enew", self.Eold, self.Enew
        return self.metropolis.acceptReject(self.Eold, self.Enew)
    
    def checkAccepted(self, quenchedE, quenched_coords, accepted):
        """
        if the step was accepted, update self.Eold
        """
        if accepted:
            self.Eold = self.Enew


