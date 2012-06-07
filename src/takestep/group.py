'''
Created on Jun 7, 2012

@author: vr274
'''

class GroupSteps(object):
    def __init__(self, steptakers):
        self.steptakers = steptakers

    def takeStep(self, coords):
        for step in self.steptakers:
            step.takeStep(coords)
       
    def updateStep(self, accepted):
        for step in self.steptakers:
            step.updateStep(accepted)

class BlockMoves(object):
    def __init__(self):
        self._steptakers = []
        self._current = 0
        self._counter = 0

    def addBlock(self, nsteps, takestep):
        self._steptakers.append([nsteps, takestep])
      
    def takeStep(self, coords):
        self._counter += 1
        if(self._counter > self._steptakers[self._current][0]):
            self._current += 1
            if(self._current >= self.steptakers.size):
                self._current = 0
        self.steptakers[self._current].takeStep(coords)
       
    def updateStep(self, accepted):
        self.steptakers[self._current].updateStep(accepted)
    
    