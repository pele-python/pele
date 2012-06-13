import numpy as np


class EnergyHistogram(object):
    """
    this class will build 1 dimensional histogram.  
    It's designed for energies, but it could work for any float data
    """
    def __init__(self, emin, emax, nbins = 100):
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        self.de = (self.emax-self.emin)/self.nbins

        self.visits = np.zeros(self.nbins)
        self.count = 0
    
    def insert(self, e):
        if not self.emin <= e < self.emax:
            print "histogram> warning: energy out of range", e
            return
        i = int((e-self.emin)/self.de)
        self.visits[i] += 1
        self.count += 1
        
    def insertWrapper(self, e, coords, acceptstep):
        return self.insert(e)

    def __iter__(self):
        return HistIter(self)
    
    def __call__(self, e, coords, acceptstep):
        return self.insert(e)
    

class HistIter(object):
    def __init__(self, hist):
        self.hist = hist
        self.counter = -1
    
    def __iter__(self):
        return self
    
    def next(self):
        self.counter += 1
        if self.counter >= self.hist.nbins:
            raise StopIteration
        return (self.hist.emin + self.hist.de * self.counter), \
            self.hist.visits[self.counter]
