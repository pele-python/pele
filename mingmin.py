import numpy as np #to access np.exp() not built int exp
import numpy.random as RNG #to access np.exp() not built int exp
from math import *
import getopt, sys
import lj
#import blj
import scipy.optimize.lbfgsb
#import ljpshift as ljpshift
import ljpshiftfast as ljpshift
#import steepest_descent
import copy
import mykeyword
import saveit

def adjustCenterOfMass(coords, natoms):
    CoM = [0.,0.,0.]
    for i in range(natoms):
        for k in range(3):
            CoM[k] += coords[i*3+k]
    for k in range(3):
        CoM[k] /= natoms
    for i in range(natoms):
        for k in range(3):
            coords[i*3+k] -= CoM[k]

def printxyz(fout, coords, natoms, E=""):
    adjustCenterOfMass(coords, natoms)
    fout.write( str(natoms) + "\n")
    fout.write( str(E) + "\n")
    for i in range(natoms):
        fout.write( "LA "+ str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 

def printcoords(fout, coords, natoms):
    for i in range(natoms):
        fout.write( str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 


def quench(potential, coords, natoms):
    initE, initV = potential.getEnergyGradient(coords)

    #with open("oldcoords", "w") as fout:
        #printxyz(fout, coords, natoms)

    #newcoords, newE = steepest_descent.steepestDescent(potential.getEnergyGradient, coords, 100)
    newcoords, newE, dictionary = scipy.optimize.fmin_l_bfgs_b(potential.getEnergyGradient, coords, iprint=-1, pgtol=1e-3)

    #with open("oldcoords", "a") as fout:
        #printxyz(fout, newcoords, natoms)

    V = dictionary["grad"]
    funcalls = dictionary["funcalls"]
    warnflag = dictionary['warnflag']
    if warnflag > 0:
        print "warning: problem with quench: ",
        if warnflag == 1:
            print "too many function evaluations"
        else:
            print dictionary['task']

    #print "quench: Ei", initE, " Ef", newE, "max Vi ", np.max(initV), "max Vf", np.max(V)
    print "quench: Ef=", newE, "steps=", funcalls, "max(V)=", np.max(np.abs(V))
    return newcoords, newE


def mcStep(potential, coordsold, natoms, Equench_old, temperature, stepsize):
    """take one monte carlo basin hopping step"""
    #########################################################################
    #take step
    #########################################################################
    coords = copy.copy(coordsold) #make  a working copy
    for j in range(natoms*3):
        rand = 2.*RNG.rand()-1.
        #print "rand ", rand
        coords[j] += stepsize*rand

    #########################################################################
    #quench
    #########################################################################
    qcoords, Equench = quench (potential, coords, natoms)

    #########################################################################
    #check whether step is accepted using metropolis algorithm.
    #Use quenched energies.
    #########################################################################
    acceptstep = True
    wcomp = (Equench - Equench_old)/temperature
    w=min(1.0,exp(-wcomp))
    rand = RNG.rand()
    if (rand > w): acceptstep = False

    print "mc step: Eo", Equench_old, "Ef", Equench, "accepted", acceptstep

    if acceptstep:
        return acceptstep, qcoords, Equench
    else:
        return acceptstep, coordsold, Equench_old

class manageStepSize:
    """a class to manage the adaptive step size"""
    def __init__(self, stepsize, accrat, nstepsaccrat, f):
        self.stepsize = stepsize
        self.f = f
        self.accrat = accrat
        self.nstepsaccrat = nstepsaccrat

        self.nsteps = 0
        self.nstepstot = 0
        self.naccepted = 0
        self.nadjust = 0
        self.changehist = []

    def getStepSize(self):
        return self.stepsize

    def adjustFactor(self, factor):
        """this is the function which ensures f is a reasonable number"""
        changef = 1.5
        if len(self.changehist) == 0: return 
        if (factor < 1.0) != (self.changehist[-1] < 1.0):
            self.f -= (self.f - 1)/1.5


    def adjustStep(self):
        """adjust stepsize"""
        self.nadjust += 1
        rat = float(self.naccepted)/self.nsteps
        if rat < self.accrat:
            #reduce step size
            self.adjustFactor( 1./self.f )
            self.stepsize /= self.f
            self.changehist.append( 1./self.f )
        else:
            #increase step size
            self.adjustFactor( self.f )
            self.stepsize *= self.f
            self.changehist.append( self.f )
        self.nsteps = 0
        self.nsteps = 0
        self.naccepted = 0
        print "accrat was ", rat, "new stepsize is ", self.stepsize, "f is", self.f

    def insertStep(self, accepted ):
        """tell us whether a step was accepted or rejected"""
        if accepted: self.naccepted += 1
        self.nsteps += 1
        self.nstepstot += 1
        if self.nsteps == self.nstepsaccrat:
            self.adjustStep()


def monteCarlo(potential, coords, natoms, nsteps, temperature, stepsize, nstepsequil, savelowest):
    fout = open("dump.q.xyz", "w")
    #########################################################################
    #do initial quench
    #########################################################################
    print "starting monteCarlo, nsteps", nsteps
    print "calculating initial energy"
    potel = potential.getEnergy(coords)
    savelowest.insert(potel, coords)
    printxyz(fout, coords, natoms, potel)
    print "initial energy", potel
    potel, V = potential.getEnergyGradient(coords)
    print "max gradient", np.max(V), potel
    print "minimizing initial coords"

    #exit()
    newcoords, newE = quench(potential, coords, natoms)
    Equench, V = potential.getEnergyGradient(newcoords)
    print "newcoords max V", np.max(V), Equench

    coords = newcoords
    printxyz(fout, coords, natoms, Equench)

    savelowest.insert(Equench, coords)
    
    #########################################################################
    #do equilibration run.  Adjust the stepsize every nstepaccrat to ensure the
    #acceptance ratio is met
    #########################################################################
    faccrat = 1.5 #the ratio by which to multiply or dived the stepsize
    accrat = 0.5
    nstepsaccrat = 10
    manstep = manageStepSize (stepsize, accrat, nstepsaccrat, faccrat)
    for istep in range(nstepsequil):
        print "step number ", istep
        acceptstep, newcoords, Equench_new = mcStep(potential, coords, natoms, Equench, temperature, manstep.stepsize)
        manstep.insertStep(acceptstep)
        if acceptstep:
            printxyz(fout, newcoords, natoms, Equench_new)
            savelowest.insert(Equench_new, newcoords)
        coords = newcoords
        Equench = Equench_new

    #########################################################################
    #loop through monte carlo steps
    #########################################################################
    stepsize = manstep.stepsize
    for istep in range(nsteps):
        print "step number ", istep + nstepsequil
        acceptstep, newcoords, Equench_new = mcStep(potential, coords, natoms, Equench, temperature, stepsize)
        if acceptstep:
            printxyz(fout, newcoords, natoms, Equench_new)
            savelowest.insert(Equench_new, newcoords)
        coords = newcoords
        Equench = Equench_new


    fout.close()



def main():
    #########################################################################
    #load coords, and get natoms
    #########################################################################
    coords2d = np.genfromtxt('coords')
    natoms = len(coords2d[:,0])
    coords = np.zeros(3*natoms, np.float64)
    for i in range(natoms):
        for k in range(3):
            coords[i*3+k] = coords2d[i,k]


    #########################################################################
    #load keywords and set defaults
    #########################################################################
    keys = mykeyword.myKeywordClass()
    keys.setDefaults()
    with open("data", "r") as fin:
        keys.readKeywords(fin)
    nstepsequil = 100

    #########################################################################
    #initialize potential, etc
    #########################################################################
    if keys.potential == "lj":
        potential = lj.LJ(keys.eps, keys.sig, natoms, keys.boxl)
    if keys.potential == "ljcpp":
        import ljcpp
        potential = ljcpp.LJ()
    if keys.potential == "binary":
        potential = ljpshift.LJpshift( natoms, keys.ntypeA, keys.boxl, keys.cutoff, keys.epsBB, keys.sigBB, keys.epsAB, keys.sigAB)
        #potential = ljpshift.LJpshift( natoms, keys.ntypeA, keys.boxl, keys.cutoff, keys.epsBB, keys.sigBB, keys.epsAB, keys.sigAB)

    #########################################################################
    #run monte carlo
    #########################################################################
    savelowest = saveit.saveit()
    monteCarlo(potential, coords, natoms, keys.nmcsteps, keys.temperature, keys.stepsize, nstepsequil, savelowest)

    #########################################################################
    #print results
    #########################################################################
    with open("lowest", "w") as fout:
        printxyz(fout,  savelowest.lowestcoords, natoms, savelowest.lowestE)

if __name__ == "__main__":
    main()
