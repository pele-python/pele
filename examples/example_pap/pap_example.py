'''
Created on Apr 16, 2012

@author: vr274
'''

'''
Created on 13 Apr 2012

@author: ruehle
'''

import potentials.potential
import pygmin
import numpy as np
import basinhopping as bh
import take_step.random_displacement as ts

class PatchyParticle(potentials.potential.potential):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
    
    def getEnergy(self, coords):
        x = self.toReal(coords).copy()
        #x = coords
        #print x.reshape(17,3)
        return pygmin.getPAPEnergy(x)
    
    def getEnergyGradient(self, coords):
        x = self.toReal(coords).copy()
        xo = coords.copy()      
        #x = coords.copy() #self.toReal(coords).copy()
        grad = np.zeros(x.shape)
        E = pygmin.getPAPEnergyGradient(x, grad)
        
        m = self.getLatticeMatrix(coords)        
        natoms = (coords.size - 3)/6
        grad[0:3*natoms] = np.dot(grad[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)
        
#        eps=1e-6
#        for i in range(6*natoms,6*natoms+3):
#            xo[i]+=eps
#            grad[i]=self.getEnergy(xo)      
#            xo[i]-=2.*eps
#            grad[i]-=self.getEnergy(xo)
#            xo[i]+=eps
#            grad[i]/=2.*eps
#            grad[i]=0.
        return E, grad
    
    def getLatticeMatrix(self, coords):
        m = np.zeros([3,3])
        m[0][0] = coords[-3]       
        m[1][1] = coords[-2]
        m[2][2] = coords[-1]
        return m

    def getInverseLatticeMatrix(self, coords):
        m = np.zeros([3,3])
        m[0][0] = 1./coords[-3]       
        m[1][1] = 1./coords[-2]
        m[2][2] = 1./coords[-1]
        return m
        
    def setLatticeMatrix(self, coords, m):
        coords[-3] = m[0][0]
        coords[-2] = m[1][1]
        coords[-1] = m[2][2]
        return m
        
    def toReduced(self, coords):
        natoms = (coords.size - 3)/6
        m = self.getInverseLatticeMatrix(coords)
        x = coords.copy()
        x[0:3*natoms] = np.dot(coords[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)
        x[-3:-1] = coords[-3:-1]
        return x
    
    def toReal(self, coords):
        natoms = (coords.size - 3)/6
        m = self.getLatticeMatrix(coords)
        x = coords.copy()
        x[0:3*natoms] = np.dot(coords[0:3*natoms].reshape([natoms,3]), m).reshape(3*natoms)
        x[-3:-1] = coords[-3:-1]
        return x

class PatchyParticleFixed(potentials.potential.potential):
    '''
    classdocs
    '''
    
    def getEnergy(self, coords):
        return pygmin.getEnergy(coords)
    
    def getEnergyGradient(self, coords):
        grad = np.zeros(coords.shape)
        E = pygmin.getEnergyGradient(coords, grad)
        return E, grad

        
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
        
def run_tests():
    pygmin.init()
    tmp = np.loadtxt('coords')
    
    x = tmp.reshape(tmp.size)
    #print tmp[1:8,:]
    pot = PatchyParticle()    
    #x[-3] -= 1e-6
#    x[0:8*6] = x[0:8*6] np.random.random(x[0:8*6].shape)*0.01
    #x[0]+=0.1
    x =pot.toReduced(x) #+ np.random.random(x.shape)*0.01
    #x[-2] += 0.5
    #x[-3] += 100
    a = x.copy()
    e=pot.getEnergy(x)
    e,g = pot.getEnergyGradient(x)
    gn = pot.NumericalDerivative(a, 1e-6)
    print e
    print g
    print gn
    #exit()
    #x[-1] += 0.1
    e,g = pot.getEnergyGradient(x)
    print "start steep"
    #steep(pot, x)
    #exit()
    
    gn = pot.NumericalDerivative(x, 1e-6)
    step = ts.takeStep( stepsize=0.1)
    opt = bh.BasinHopping(x, pot, takeStep=step.takeStep, quenchRoutine = steep)
    opt.run(10)
    #print pot.getEnergyGradient(opt.coords)
    
    

if __name__ == '__main__':
    run_tests()