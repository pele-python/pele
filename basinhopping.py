# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import numpy.random as RNG


class BasinHopping:
  def __init__(self, coords, potential, takestep, temperature=1.0, storage=None, manstep=None):
    self.coords = coords
    self.temperature = temperature
    self.storage = storage
    self.potential = potential
    self.takestep = takestep
    self.manstep = manstep

  def run(self, nsteps):
    #########################################################################
    #do initial quench
    #########################################################################
    print "doing monteCarlo, nsteps", nsteps
    potel = self.potential.getEnergy(self.coords)
    print "initial energy", potel

    if(self.storage):
      self.storage(potel, self.coords)
      
    potel, V = self.potential.getEnergyGradient(self.coords)
    print "max gradient", np.max(V), potel
    print "minimizing initial coords"

    newcoords, newE = self.quench(self.coords)
    Equench, V = self.potential.getEnergyGradient(newcoords)
    print "newcoords max V", np.max(V), Equench

    coords = newcoords

    if(self.storage):
      self.storage(Equench, coords)
    
    for istep in xrange(nsteps):
        print istep
        acceptstep, newcoords, Equench_new = self.mcStep(self.coords, Equench)
        if self.manstep:
          self.manstep.insertStep(acceptstep)
        if acceptstep:
            if(self.storage):
              self.storage(Equench_new, newcoords)
        self.coords = newcoords
        Equench = Equench_new
    
  def quench(self, coords):
    import scipy.optimize.lbfgsb
    #newcoords, newE = steepest_descent.steepestDescent(potential.getEnergyGradient, coords, 100)
    newcoords, newE, dictionary = scipy.optimize.fmin_l_bfgs_b(self.potential.getEnergyGradient, coords, iprint=-1, pgtol=1e-3)

    V = dictionary["grad"]
    funcalls = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    if warnflag > 0:
        print "warning: problem with quench: ",
        if warnflag == 1:
            print "too many function evaluations"
        else:
            print dictionary['task']

    print "quench: Ef=", newE, "steps=", funcalls, "max(V)=", np.max(np.abs(V))
    return newcoords, newE
    
  def mcStep(self, coordsold, Equench_old):
    """take one monte carlo basin hopping step"""
    #########################################################################
    #take step
    #########################################################################
    coords = coordsold.copy() #make  a working copy
    self.takestep.takeStep(coords)

    #########################################################################
    #quench
    #########################################################################
    qcoords, Equench = self.quench (coords)

    #########################################################################
    #check whether step is accepted using metropolis algorithm.
    #Use quenched energies.
    #########################################################################
    acceptstep = True
    wcomp = (Equench - Equench_old)/self.temperature
    w=min(1.0,exp(-wcomp))
    rand = RNG.rand()
    if (rand > w): acceptstep = False

    print "mc step: Eo", Equench_old, "Ef", Equench, "accepted", acceptstep

    if acceptstep:
        return acceptstep, qcoords, Equench
    else:
        return acceptstep, coordsold, Equench_old

