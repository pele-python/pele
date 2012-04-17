import copy

class saveit:
    """this class will keep track of the minima with the lowest energy"""
    def __init__(self, E=1000000., coords=[]):
        self.lowestE = E
        self.lowestcoords = copy.copy(coords)

    def insert(self, E, coords):
        if E < self.lowestE:
            self.lowestE = E
            self.lowestcoords = copy.copy(coords)
