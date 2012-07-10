'''
Created on Jun 7, 2012

@author: vr274
'''

class GroupSteps(object):
    def __init__(self, steptakers):
        self.steptakers = steptakers

    def takeStep(self, coords, **kwargs):
        for step in self.steptakers:
            step.takeStep(coords, **kwargs)
       
    def updateStep(self, accepted, **kwargs):
        for step in self.steptakers:
            step.updateStep(accepted, **kwargs)

class BlockMoves(object):
    def __init__(self):
        self._steptakers = []
        self._current = 0
        self._counter = 0

    def addBlock(self, nsteps, takestep):
        self._steptakers.append([nsteps, takestep])
      
    def takeStep(self, coords, **kwargs):
        self._counter += 1
        if(self._counter > self._steptakers[self._current][0]):
            self._current += 1
            self._counter = 0
            if(self._current >= len(self._steptakers)):
                self._current = 0
        self._steptakers[self._current][1].takeStep(coords, **kwargs)
       
    def updateStep(self, accepted, **kwargs):
        self._steptakers[self._current][1].updateStep(accepted, **kwargs)

class Reseeding(object):
    def __init__(self, takestep, reseed, maxnoimprove=100):
        self.takestep = takestep
        self.reseed = reseed
        self.maxnoimprove = maxnoimprove
        self._noimprove = 0
        self.lowest = None
        
    def takeStep(self, coords, **kwargs):
        if self._noimprove >= self.maxnoimprove:
            print "The energy did not improve after " + str(self._noimprove) + \
                " steps, reseeding"
            self.reseed.takeStep(coords, **kwargs)
            kwargs['driver'].acceptTest.forceAccept()
        else:
            self.takestep.takeStep(coords, **kwargs)
       
    def updateStep(self, accepted, **kwargs):
        driver = kwargs["driver"]
        if self.lowest == None:
            self.lowest = driver.markovE
        if self._noimprove >= self.maxnoimprove:
            self.reseed.updateStep(accepted, **kwargs)
            self._noimprove=1
        else:
            self.takestep.updateStep(accepted, **kwargs)
            if driver.markovE >= self.lowest:
                self._noimprove+=1
            else:
                self.lowest = driver.markovE
                self._noimprove=1
            