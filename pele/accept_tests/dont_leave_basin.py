import pele.exceptions as exc

__all__ = ["DontLeaveBasin"]


class DontLeaveBasin(object):
    """
    reject the step if the new energy is different from the old energy
    """

    def __init__(self, Ecriterion=1e-4):
        if Ecriterion < 0:
            raise exc.SignError("Energy criterion must be positive.")
        self.Ecriterion = Ecriterion

    def acceptReject(self, Eold, Enew, **kwargs):
        return abs(Eold - Enew) < self.Ecriterion
