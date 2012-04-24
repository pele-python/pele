import numpy as np
import copy
import rotations as rot
import potentials.potential as potential

class MinPermDistPotential(potential.potential):
    def __init__(self, XA, XB, L=0.4):
        self.XA0 = copy.copy(XA)
        self.XB0 = copy.copy(XB)
        self.XB = copy.copy(XB)
        self.nsites = len(XA)/3
        self.L2 = L*L
        pass

    def getEnergy(self, AA ):
        # Rotate XB0 according to angle axis AA
        rot_mx = rot.aa2mx( AA )
        for j in range(self.nsites):
            i = 3*j
            self.XB[i:i+3] = np.dot( rot_mx, self.XB0[i:i+3] )
        #return distance between XB and XA0
        E = 0.
        for i1 in range(self.nsites):
            i = i1*3
            for j1 in range(self.nsites):
                j = j1*3
                r2 = np.sum( (self.XB[i:i+3] - self.XA0[j:j+3])**2 )
                E -= np.exp(-r2/self.L2)
        return E


def take_step( aa):
    aanew = rot.random_aa()
    #print "aanew ", aanew
    aa[:] = aanew[:]
    #print "aa    ", aa
