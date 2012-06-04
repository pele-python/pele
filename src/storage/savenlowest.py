'''
Created on Apr 18, 2012

@author: vr274
'''

import threading

class SaveN(object):
    '''
    Stores only the nsave lowest minima. Minima are considered as different
    if energy differs by more than accuracy
    '''

    def __init__(self, nsave=1, accuracy=1e-3, onMinimumAdded=None, onMinimumRemoved=None):
        '''
        Constructor
        '''
        self.nsave=nsave
        self.data=[]
        self.accuracy=accuracy
        self.onMinimumAdded=onMinimumAdded
        self.onMinimumRemoved=onMinimumRemoved
        self.next_free_id = 1
        self.lock = threading.Lock()
        
    def insert(self, E, coords):
        # does minima already exist, if yes exit?
        self.lock.acquire()
        for i in self.data:
            if(abs(i[0] - E) < self.accuracy):
                self.lock.release()
                return
                
        id = self.next_free_id
        self.next_free_id+=1
        # otherwise, add it to list & sort
        new = (E, coords.copy(), id)
        self.data.append(new)
        self.data.sort()
        if(self.onMinimumAdded):
            self.onMinimumAdded(new)
        # remove if too many entries
        while(len(self.data) > self.nsave):
            removed = self.data.pop()
            if(self.onMinimumRemoved):
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
        ddict = self.__dict__.copy();
        ddict["onMinimumAdded"]=None
        ddict["onMinimumRemoved"]=None
        ddict["lock"]=None
        return ddict #.items()
        
if __name__ == "__main__":
    import numpy as np
    save = SaveN(nsave = 2)
    save.insert(1., np.random.random(10))
    save.insert(3., np.random.random(10))
    save.insert(2., np.random.random(10))
    save.insert(1.001, np.random.random(10))
    for i in save.data:
        print i
        