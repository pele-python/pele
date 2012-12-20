'''
Created on 6 Apr 2012

@author: ruehle
'''

import numpy as np
import copy

__all__ = ["QuenchBenchmark"]

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
            print "Testing Minimizer " + minimizer[0]
            
            E,grad = self.potential.getEnergyGradient(coords)
            #self.potential.energies.append(E)
            x, E, tmp, tmp = minimizer[1](coords, self.potential.getEnergyGradient)
            minimizer[2] = E
            minimizer[3] = np.array(self.potential.energies).copy()-Emin
            print "Minimizer " + minimizer[0] + ": " + str(E)
            
    def plot(self):
        import pylab as pl
        for m in self.minimizer:
            pl.loglog(np.array(m[3]), label=m[0])
        
        pl.legend()
        #pl.semilogy()
        pl.xlabel("energy evaluations")
        pl.ylabel("energy")
        pl.show()

if __name__ == "__main__":
    import pygmin.potentials.lj as lj
    import scipy.optimize
    from pygmin.optimize import _quench as quench
    print "Running benchmark with lennard jones potential"
    pot = lj.LJ()
    
    natoms = 36
    
    coords = np.random.random(3*natoms)*10.0
    coords, E, tmp, tmp=quench.lbfgs_scipy(coords, pot.getEnergyGradient, tol=1e-3)
    
    coords = coords + np.random.random(coords.shape)*0.1
    tmp, Emin, tmp, tmp = quench.lbfgs_scipy(coords, pot.getEnergyGradient, tol=1e-3)
    
    bench = QuenchBenchmark(pot)
    bench.addMinimizer("lbfgs", quench.lbfgs_scipy)
    bench.addMinimizer("mylbfgs", quench.mylbfgs)
    bench.addMinimizer("lbfgs_py", quench.lbfgs_py)
    bench.addMinimizer("lbfgs_ase", quench.lbfgs_ase)
    bench.addMinimizer("cg", quench.cg)
    bench.addMinimizer("fire", quench.fire)
    bench.addMinimizer("bfgs", quench.bfgs)
    #bench.addMinimizer("fmin", quench.fmin)
    #bench.addMinimizer("steep", quench.steepest_descent)
    
    print "The reference energy is " + str(Emin)
    bench.run(Emin,coords)
    bench.plot()
