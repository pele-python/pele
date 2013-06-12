import numpy as np
from generic import Takestep

class RandomRotation(Takestep):
    def __init__(self, stepsize=1.0):
        Takestep.__init__(self, stepsize=stepsize)
    def takeStep(self, coords, **kwargs):
        rotation_matrix = np.random.uniform(low=-self.stepsize,
                                      high=self.stepsize, 
                                      size=coords[self.srange].shape)