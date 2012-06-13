'''
Created on 6 Apr 2012

@author: ruehle
'''

import numpy as np
import copy

class PotentialWrapper(object):
    def __init__(self, potential):
        self.potential=potential
        self.reset()
        
    def reset(self):
        self.energies=[]
        
    def getEnergy(self, coords):
        E = self.potential.getEnergy(coords)
        self.energies.append(E)
        return E    

    def getEnergyGradient(self, coords):
        E,g = self.potential.getEnergyGradient(coords)
        self.energies.append(E)
        return E,g
     
    def getGradient(self, coords):
        E,g = self.potential.getEnergyGradient(coords)
        self.energies.append(E)
        return g
        
class QuenchBenchmark(object):
    '''
    classdocs
    '''


    def __init__(self, potential):
        '''
        Constructor
        '''
        self.potential=PotentialWrapper(potential)
        self.minimizer=[]
        
    def addMinimizer(self, label, minimizer):
        self.minimizer.append([label, minimizer, 0.0, None])
        
    def run(self, Emin,coords):
        for minimizer in self.minimizer:
            self.potential.reset()
            E,x = minimizer[1](self.potential, coords)
            minimizer[2] = E
            minimizer[3] = np.array(self.potential.energies).copy()-Emin
            print E,Emin
            
    def plot(self):
        import pylab as pl
        for m in self.minimizer:
            pl.loglog(np.array(m[3]), label=m[0])
        
        pl.legend()
        #pl.semilogy()
        pl.xlabel("energy evaluations")
        pl.ylabel("energy")
        pl.show()
    
    
def lbfgs(pot, coords):
    tmp,E,dictio = scipy.optimize.fmin_l_bfgs_b(pot.getEnergyGradient, coords, iprint=-1, pgtol=1e-8)
    return E,tmp

def lbfgs2(pot, coords):
    tmp,E,dictio = scipy.optimize.fmin_l_bfgs_b(pot.getEnergy, coords, fprime=pot.getGradient, iprint=-1, pgtol=1e-6)
    return E,tmp

def cg(pot, coords):
    r = scipy.optimize.fmin_cg(pot.getEnergy, np.array(coords), fprime=pot.getGradient)
    return r[1],r[0]

if __name__ == "__main__":
    import potentials.lj as lj
    import scipy.optimize
    
    print "Running benchmark with lennard jones potential"
    pot = lj.LJ()
    
    natoms = 36
    
    E,coords=lbfgs(pot, np.random.random(3*natoms)*10.0)
    coords = coords + np.random.random(coords.shape)*0.2
    Emin,tmp = lbfgs(pot,coords)
    
    bench = QuenchBenchmark(pot)
    bench.addMinimizer("lbfgs", lbfgs2)
    bench.addMinimizer("cg", cg)
    print Emin
    bench.run(Emin,coords)
    bench.plot()
