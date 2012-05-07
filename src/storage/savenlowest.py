'''
Created on Apr 18, 2012

@author: vr274
'''

class SaveN(object):
    '''
    Stores only the nsave lowest minima. Minima are considered as different
    if energy differs by more than accuracy
    '''

    def __init__(self, nsave=1, accuracy=1e-3, onNewMinimumFound=None):
        '''
        Constructor
        '''
        self.nsave=nsave
        self.data=[]
        self.accuracy=accuracy
        self.onNewMinimumFound=onNewMinimumFound
    
    def insert(self, E, coords):
        # does minima already exist, if yes exit?
        for i in self.data:
            if(abs(i[0] - E) < self.accuracy):
                return
        
        # otherwise, add it to list & sort
        self.data.append((E, coords.copy()))
        icoords = self.data[-1][1]
        self.data.sort()
        # remove if too many entries
        while(len(self.data) > self.nsave):
            self.data.pop()
        if(self.onNewMinimumFound):
            self.onNewMinimumFound(E, icoords)
            
if __name__ == "__main__":
    import numpy as np
    save = SaveN(nsave = 2)
    save.insert(1., np.random.random(10))
    save.insert(3., np.random.random(10))
    save.insert(2., np.random.random(10))
    save.insert(1.001, np.random.random(10))
    for i in save.data:
        print i
        