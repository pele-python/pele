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

def adjustCenterOfMass(coords, natoms):
    CoM = np.zeros(3, np.float64)
    for i in range(natoms):
        CoM += coords[i*3:i*3+3]
    CoM /= natoms
    for i in range(natoms):
        coords[i*3:i*3+3] -= CoM

def printxyz(fout, coords, natoms, E=""):
    adjustCenterOfMass(coords, natoms)
    fout.write( str(natoms) + "\n")
    fout.write( str(E) + "\n")
    for i in xrange(natoms):
        fout.write( "LA "+ str(coords[i*3+0])+" "+ str(coords[i*3+1])+" "+ str(coords[i*3+2])+" "+ "\n" ) 

def printcoords(fout, coords, natoms):
    for i in xrange(natoms):
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


def mcStep(potential, coordsold, natoms, Equench_old, temperature, takeStep):
    """take one monte carlo basin hopping step"""
    #########################################################################
    #take step
    #########################################################################
    coords = copy.copy(coordsold) #make  a working copy
    takeStep.takeStep(coords)

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


def monteCarlo(potential, coords, natoms, nsteps, temperature, stepsize, nstepsequil, savelowest, accrat=0.5, accrat_frq=50):
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
    manstep = adaptive_step.manageStepSize (stepsize, accrat, accrat_frq)
    takeStep = take_step.takeStep( RNG = np.random.rand, getStep = manstep.getStepSize )
    for istep in xrange(nstepsequil):
        print "step number ", istep
        acceptstep, newcoords, Equench_new = mcStep(potential, coords, natoms, Equench, temperature, takeStep )
        manstep.insertStep(acceptstep)
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
    #run monte carlo
    #########################################################################
    savelowest = saveit.saveit()
    monteCarlo(potential, coords, natoms, keys.nmcsteps, keys.temperature, keys.stepsize, nstepsequil, savelowest, keys.accrat, keys.accrat_frq )

    #########################################################################
    #print results
    #########################################################################
    with open("lowest", "w") as fout:
        printxyz(fout,  savelowest.lowestcoords, natoms, savelowest.lowestE)

if __name__ == "__main__":
    main()
