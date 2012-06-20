class TakestepInterface(object):
    def takeStep(self, coords, **kwargs):
        pass
        
    def updateStep(self, accepted, **kwargs):
        pass
    
    def scale(self, factor):
        pass
    
class Takestep(TakestepInterface):
    def __init__(self, stepsize=1.0):
        self.stepsize = stepsize
        
    def scale(self, factor):
        self.stepsize*=factor
        
class TakestepSlice(Takestep):
    def __init__(self, srange=None, stepsize=1.0):
        self.stepsize=stepsize
        self.srange = srange