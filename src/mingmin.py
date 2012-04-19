import numpy as np #to access np.exp() not built int exp
from math import *
import getopt, sys
import scipy.optimize.lbfgsb
import copy
import mykeyword
import basinhopping as bh
import printing.print_atoms_xyz as printxyz

def adjustCenterOfMass(coords):
    natoms = len(coords)/3
    CoM = np.zeros(3, np.float64)
    for i in range(natoms):
        CoM += coords[i*3:i*3+3]
    CoM /= natoms
    for i in range(natoms):
        coords[i*3:i*3+3] -= CoM

class printWrapper:
    def __init__(self, fout):
        self.fout = fout
    def printIt(self, E, coords, accept):
        if accept:
            adjustCenterOfMass(coords)
            printxyz.printAtomsXYZ(self.fout, coords, E)

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



    #########################################################################
    #initialize basin hopping class
    #########################################################################

    #now initialize the basin hopping class
    opt = keys.getBasinHopping( coords, natoms )

    #optional ...
    #class to dump xyz files after each successful step
    #add to list of event_after_step
    fdump = open("dump.q.xyz", "w")
    dumper = printWrapper(fdump)
    opt.event_after_step.append(dumper.printIt)

    #now run basin hopping
    opt.run(keys.nmcsteps)
    fdump.close()

    #########################################################################
    #print results
    #########################################################################
    print "Finished:"
    with open("lowest", "w") as fout:
        for (E,X) in keys.savelowest.data:
            adjustCenterOfMass(X)
            printxyz.printAtomsXYZ(fout,  X, E)

    print "Acceptance ratio for run= ", keys.manstep.accratTot(), "Step=", keys.takeStep.getStep()
    print "lowest energy found ", keys.savelowest.data[0][0]

if __name__ == "__main__":
    main()
