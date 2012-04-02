# -*- coding: iso-8859-1 -*-
import numpy as np
import scipy
from math import *
import metropolis


class BasinHopping:
  """A class to run the basin hopping algorithm
  
     input:

     coords: the initial set of coordinates

     takestep: the function which randomly perterbs the system, e.g. random
                dispacement.  It takes the form

            takestep(coords)

     acceptTests: a list of functions which return false if a quench should be rejected
                   Each test in the list takes the form

            test(Eold, Enew, new_coords):

     event_after_step: an optional list of functions which act just after each monte carlo
                        round.  Each even in the list takes the form

            event(Equench_new, newcoords, acceptstep)
  """
  def __init__(self, coords, potential, takestep, storage=None, manstep=None, event_after_step=[], \
          acceptTests=[]  \
          ):
    self.coords = coords
    self.storage = storage
    self.potential = potential
    self.takestep = takestep
    self.manstep = manstep
    self.event_after_step = event_after_step
    self.acceptTests = acceptTests
    if len(self.acceptTests) == 0: #set default value
        self.metrop_test = metropolis.Metropolis(1.0)
        self.acceptTests=[self.metrop_test.acceptReject]

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
        for event in self.event_after_step:
            event(Equench_new, newcoords, acceptstep)
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
    #check whether step is accepted with user defined tests.  If any returns
    #false then reject step.
    #########################################################################
    acceptstep = True
    for test in self.acceptTests:
        if not test(Equench_old, Equench, qcoords):
            acceptstep = False

    #########################################################################
    #return
    #########################################################################
    if acceptstep:
        return acceptstep, qcoords, Equench
    else:
        return acceptstep, coordsold, Equench_old

