__all__ = ["DontLeaveBasin"]

class DontLeaveBasin:
    """
    reject the step if the new energy is different from the old energy
    """
    def __init__(self, Ecriterion = 1e-4):
        if Ecriterion < 0:
            raise SignError, "Energy criterion must be positive."
        self.Ecriterion = Ecriterion

    def acceptReject(self, Eold, Enew, qcoords=[], coords=[]):
        return abs(Eold - Enew) < self.Ecriterion
    
class SignError(Exception):
    """
    this is an exception to be raised if the energy criterion is negative
    """
    pass