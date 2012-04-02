import numpy as np #to access np.exp() not built int exp
import numpy.random as RNG #to access np.exp() not built int exp
from math import *
import getopt, sys
import scipy.optimize.lbfgsb
import copy
import mykeyword
import saveit
import adaptive_step
import take_step
import basinhopping as bh
import metropolis

def adjustCenterOfMass(coords, natoms):
    CoM = np.zeros(3, np.float64)
    for i in range(natoms):
        CoM += coords[i*3:i*3+3]
    CoM /= natoms
    for i in range(natoms):
        coords[i*3:i*3+3] -= CoM

def printxyz(fout, coords, E=""):
    natoms = len(coords)/3
    adjustCenterOfMass(coords, natoms)
    fout.write( str(natoms) + "\n")
    fout.write( str(E) + "\n")
    for i in xrange(natoms):
        fout.write( "LA "+ str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 

def printcoords(fout, coords, natoms):
    for i in xrange(natoms):
        fout.write( str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 

class printWrapper:
    def __init__(self, fout):
        self.fout = fout
    def printIt(self, E, coords, accept):
        if accept:
            printxyz(self.fout, coords, E)

def main():
    #########################################################################
    #load coords, and get natoms
    #########################################################################
    coords2d = np.genfromtxt('coords')
    natoms = len(coords2d[:,0])
    coords = np.zeros(3*natoms, np.float64)
    for i in xrange(natoms):
        coords[i*3:i*3+3] = coords2d[i,0:3]


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
        import potentials.lj as lj
        potential = lj.LJ(keys.eps, keys.sig, keys.boxl)
    if keys.potential == "ljcpp":
        import potentials.ljcpp as ljcpp
        potential = ljcpp.LJ()
    if keys.potential == "binary":
        usefortran = True
        if usefortran:
            import potentials.ljpshiftfast as ljpshift
        else:
            import potentials.ljpshift as ljpshift
        potential = ljpshift.LJpshift( natoms, keys.ntypeA, keys.boxl, keys.cutoff, keys.epsBB, keys.sigBB, keys.epsAB, keys.sigAB)

    #########################################################################
    #initialize basing hopping class
    #########################################################################
    savelowest = saveit.saveit() #class to save the lowest energy structure
    manstep = adaptive_step.manageStepSize (keys.stepsize, keys.accrat, keys.accrat_frq) #class to do step size adjustment
    takeStep = take_step.takeStep( RNG = np.random.rand, getStep = manstep.getStepSize ) #class to impliment the take step routine
    #class to dump xyz files after each successful step
    fdump = open("dump.q.xyz", "w")
    dumper = printWrapper(fdump)
    event_after_step = [dumper.printIt]
    #class to impiment acceptence criterion
    metrop_test = metropolis.Metropolis(keys.temperature)
    acceptTests=[metrop_test.acceptReject]
    #now initialize the basin hopping class
    opt = bh.BasinHopping(coords, potential, takeStep, storage = savelowest.insert, manstep = manstep, \
            event_after_step=event_after_step, \
            acceptTests=acceptTests, \
            )
    opt.run(keys.nmcsteps)
    fdump.close()

    #########################################################################
    #print results
    #########################################################################
    with open("lowest", "w") as fout:
        printxyz(fout,  savelowest.lowestcoords, savelowest.lowestE)

if __name__ == "__main__":
    main()
