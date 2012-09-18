class DontLeaveBasin:
    """
    reject the step if the new energy is different from the old energy
    """
    def __init__(self, Ecriterion = 1e-4):
        self.Ecriterion = Ecriterion

    def acceptReject(self, Eold, Enew, qcoords=[], coords=[]):
        return abs(Eold - Enew) < self.Ecriterion
