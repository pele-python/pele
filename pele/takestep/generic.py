__all__ = ["TakestepInterface", "Takestep", "TakestepSlice"]


class TakestepInterface(object):
    """Interface for step taking classes"""

    def takeStep(self, coords, **kwargs):
        """take a step

        Parameters
        ----------
        coords : array like object
            coordinates
        """
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        return self.takeStep(*args, **kwargs)

    def updateStep(self, accepted, **kwargs):
        """feedback from basin hopping if last step was accepted

        Parameters
        ----------

        accepted : boolean
            True if the last step was accepted, otherwise false
        """
        pass

    def scale(self, factor):
        """scale the stepsize

        This function needs to be implemented to support adaptive step taking

        Parameters
        ----------
        factor : float
            factor to scale the stepsize
        """
        pass


class Takestep(TakestepInterface):
    """basic takestep interface which stores the stepsize"""

    def __init__(self, stepsize=1.0):
        self.stepsize = stepsize

    def scale(self, factor):
        self.stepsize *= factor


class TakestepSlice(Takestep):
    """basic takestep interface on slice of coordinates array"""

    def __init__(self, srange=None, stepsize=1.0):
        self.stepsize = stepsize
        self.srange = srange
