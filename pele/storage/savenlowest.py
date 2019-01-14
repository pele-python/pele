"""
Created on Apr 18, 2012

@author: vr274
"""
from __future__ import print_function

import threading
import operator

class Minimum(object):
    """
    class for storing minima
    """
    
    def __init__(self, E, coords):        
        self.energy = E
        self.coords = coords.copy()        
    
class SaveN(object):
    """
    Stores only the nsave lowest minima. Minima are considered as different
    if energy differs by more than accuracy
    """

    def __init__(self, nsave=1, accuracy=1e-3, onMinimumAdded=None, onMinimumRemoved=None, compareMinima=None):
        """
        Constructor
        """
        self.nsave=nsave
        self.data=[]
        self.accuracy=accuracy
        self.onMinimumAdded=onMinimumAdded
        self.onMinimumRemoved=onMinimumRemoved
        self.compareMinima=compareMinima
        self.lock = threading.Lock()
        
    def __call__(self, E, coords):
        self.insert(E, coords)

    def addMinimum(self, E, coords):
        self.insert(E, coords)
        
    def insert(self, E, coords):
        new = Minimum(E, coords)
        # does minima already exist, if yes exit?
        self.lock.acquire()
        for i in self.data:
            if abs(i.energy - E) < self.accuracy:
                if self.compareMinima:
                    if not self.compareMinima(new, i):
                        continue                
                
                self.lock.release()
                return
                
        # otherwise, add it to list & sort
        self.data.append(new)
        self.data.sort(key=operator.attrgetter('energy'))
        if self.onMinimumAdded:
            self.onMinimumAdded(new)
        # remove if too many entries
        while len(self.data) > self.nsave:
            removed = self.data.pop()
            if self.onMinimumRemoved:
                self.onMinimumRemoved(removed)
        self.lock.release()
        
    def save(self, filename):
        import pickle
        output = open(filename, "w")
        pickle.dump(self, output)
        
    @classmethod
    def load(filename):
        import pickle
        infile = open(filename, "w")
        return pickle.load(infile)
    
    
    def __getstate__(self):
        ddict = self.__dict__.copy()
        ddict["onMinimumAdded"]=None
        ddict["onMinimumRemoved"]=None
        del ddict["lock"]
        return ddict #.items()
    
    def __setstate__(self, dct):
        self.__dict__.update(dct)
        self.lock = threading.Lock()
        
    
        
if __name__ == "__main__":
    import numpy as np
    save = SaveN(nsave = 2)
    save.insert(1., np.random.random(10))
    save.insert(3., np.random.random(10))
    save.insert(2., np.random.random(10))
    save.insert(1.001, np.random.random(10))
    for i in save.data:
        print(i.energy, i.coords)
        
