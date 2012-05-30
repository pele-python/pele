import numpy as np


class EnergyHistogram(object):
    def __init__(self, emin, emax, nbins):
        self.emin = emin
        self.emax = emax
        self.nbins = nbins
        self.de = (self.emax-self.emin)/self.nbins
        
        self.visits = np.zeros(nbins)
        self.count = 0
    
    def insert(self, e):
        if not self.emin <= e < self.emax:
            print "histogram> warning: energy out of range", e
            return
        i = int((e-self.emin)/self.de)
        self.visits[i] += 1
        self.count += 1
        
    def insertWrapper(self, e, crap1, crap2):
        return self.insert(e)
    
    def __iter__(self):
        return HistIter(self)
    

class HistIter(object):#
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