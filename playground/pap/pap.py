import potentials.potential
import numpy as np
import pele


class PatchyParticle(potentials.potential.potential):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
    
    def getEnergy(self, coords):
        x = self.toReal(coords)
        return pele.getPAPEnergy(x)
    
    def getEnergyGradient(self, coords):
        x = self.toReal(coords)
        grad = np.zeros(x.shape)
        E = pele.getPAPEnergyGradient(x, grad)

        # transform gradient to reduced units
        m = self.getLatticeMatrix(coords)        
        natoms = (coords.size - 3)/6
        grad[0:3*natoms] = np.dot(grad[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)

# the lower part would calculate numerical lattice gradient        
#        eps=1e-6
#        for i in range(6*natoms,6*natoms+3):
#            xo[i]+=eps
#            grad[i]=self.getEnergy(xo)      
#            xo[i]-=2.*eps
#            grad[i]-=self.getEnergy(xo)
#            xo[i]+=eps
#            grad[i]/=2.*eps
        #print "bla"
        return E, grad
    
    # construct the lattice matrix
    def getLatticeMatrix(self, coords):
        m = np.zeros([3,3])
        m[0][0] = coords[-3]       
        m[1][1] = coords[-2]
        m[2][2] = coords[-1]
        return m

    # construct the inverse lattice matrix
    def getInverseLatticeMatrix(self, coords):
        m = np.zeros([3,3])
        m[0][0] = 1./coords[-3]       
        m[1][1] = 1./coords[-2]
        m[2][2] = 1./coords[-1]
        return m
        
    # set coords based on lattice matrix
    def setLatticeMatrix(self, coords, m):
        coords[-3] = m[0][0]
        coords[-2] = m[1][1]
        coords[-1] = m[2][2]
        return m
        
    # to to reduced coordinates
    def toReduced(self, coords):
        natoms = (coords.size - 3)/6
        m = self.getInverseLatticeMatrix(coords)
        x = coords.copy()
        x[0:3*natoms] = np.dot(coords[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)
        x[-3:-1] = coords[-3:-1]
        return x
    
    # go to real cartesian coordinates
    def toReal(self, coords):
        natoms = (coords.size - 3)/6
        m = self.getLatticeMatrix(coords)
        x = coords.copy()
        x[0:3*natoms] = np.dot(coords[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)
        x[-3:-1] = coords[-3:-1]
        return x


# not used, just for testing       
def steep(coords, pot):
    work = coords.copy()
    Elast = pot(work)
    step=0.0001  
    for i in xrange(10000):        
        nsteps=i
        E,grad = pot(work)
        
        dx = step*grad[:]
        if(np.max(np.abs(dx)) > 0.01):
            #step=0.01/np.max(dx)
            dx*=0.01/np.max(dx)
        #print "max", np.max(np.abs(dx))
        tmp = work-dx
        
#        while(E > Elast):
#            step*=0.5
#            print "reducing step to", step
#            tmp = work - step*grad[:]
#            tmp = work - step*grad[:]

        work = tmp
        #print E
        norm = np.linalg.norm(grad)
        if(norm < 1e-3):
             break
    #print tmp, grad
    return tmp, E, np.linalg.norm(grad), nsteps 
