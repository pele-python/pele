import numpy as np

__all__ = ["Metropolis"]


class Metropolis(object):
    """Accept steps based on the metropolis criterion
    
    Parameters
    ----------
    temperature : float
    random : callable
        return a random number in [0,1)
    """

    def __init__(self, temperature, random=np.random.rand):
        self.random = random
        self.temperature = temperature
        self._accept_next = False

    def acceptRejectE(self, Eold, Enew):
        """the Metropolis criterion"""
        if self._accept_next:
            self._accept_next = False
            return True
        if Enew < Eold: return True
        acceptstep = True
        wcomp = (Enew - Eold) / self.temperature
        w = min(1.0, np.exp(-wcomp))
        rand = self.random()
        if rand > w: acceptstep = False

        return acceptstep

    def forceAccept(self):
        """Force acceptance of the next step. This is useful for reseeding.
        """
        self._accept_next = True

    def __call__(self, Eold, Enew, qcoords, coords):
        """wrapper for acceptRejectE"""
        return self.acceptRejectE(Eold, Enew)

# class MetropolisNonQuench(object):
# """
# perform metropolis criterion on non quenched energy
# """
# def __init__(self, temperature, potential, random=np.random.rand):
# self.potential = potential
# self.metropolis = Metropolis(temperature, random)
#
#
# def acceptReject(self, Equench_old=0., Equench_new=0., qcoords=None, coords=None):
# self.Enew = self.potential.getEnergy(coords)
# try:
# self.Eold
# except AttributeError:
#            self.Eold = self.Enew
#        return self.metropolis.acceptReject(self.Eold, self.Enew)
#    
#    def checkAccepted(self, quenchedE, quenched_coords, accepted):
#        """if the step was accepted, update self.Eold
#        """
#        if accepted:
#            self.Eold = self.Enew


if __name__ == "__main__":
    met = Metropolis(1.0)
    
    
    

