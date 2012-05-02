import distpot
import numpy as np
import copy
import rotations as rot
from optimize.quench import quench
import basinhopping
import storage.savenlowest as storage
from mindistutils import CoMToOrigin, alignRotation, findBestPermutationRBMol, getDistaa


def coordsApplyRotation(coordsin, aa):
    coords = coordsin.copy()
    nmol = len(coords) / 3 / 2
    rmat = rot.aa2mx(aa)
    #rotate center of mass coords
    for imol in range(nmol):
        k = imol*3
        coords[k : k + 3] = np.dot(rmat, coords[k : k + 3])
    #update aa coords
    for imol in range(nmol):
        k = nmol*3 + imol*3
        coords[k : k + 3] = rot.rotate_aa(coords[k : k + 3], aa)
    return coords



def minPermDistRBMol(coords1, coords2, mysys, niter = 100, permlist = None):
    """
    Minimize the distance between two clusters.  The following symmetries will be accounted for
    
    Translational symmetry

    Global rotational symmetry

    Permutational symmetry

    
    This method uses basin hopping to find the rotation of X2 which best
    optimizes the overlap (an effective energy) between X1 and X2.  The overlap
    is defined to be permutation independent.

    Once the rotation is optimized, the correct permutation can be determined
    deterministically using the Hungarian algorithm.
    
    input:
    
    coords1, coords2: the structures for which to minimize the distance in com + angle-axis format
    
    permlist  ([range(natoms)]) 
        A list of lists of atoms which are interchangable.
        e.g. for a 50/50 binary mixture, permlist = [ range(1,natoms/2), range(natoms/2,natoms) ]
    """
    nmol = mysys.nmol
    
    if permlist == None:
        permlist = [range(nmol)]
    
    coords1in = coords1.copy()
    coords2in = coords2.copy()
    
    comcoords1 = copy.copy(coords1[0:3*nmol])
    aacoords1  = copy.copy(coords1[3*nmol:])
    
    comcoords2 = copy.copy(coords2[0:3*nmol])
    aacoords2  = copy.copy(coords2[3*nmol:])


    ###############################################
    # move the centers of mass to the origin
    ###############################################
    #warning: this assumes all molecules have the same mass
    coords1[0:3*nmol] = CoMToOrigin(coords1[0:3*nmol])
    coords2[0:3*nmol] = CoMToOrigin(coords2[0:3*nmol])

    coords1in = coords1.copy()
    coords2in = coords2.copy()

    comcoords1 = copy.copy(coords1[0:3*nmol])
    comcoords2 = copy.copy(coords2[0:3*nmol])


    ###############################################
    # set initial conditions
    ###############################################
    aamin = np.array([0.,0.,0.])

    ######################################
    # set up potential, dependent only on the center of mass coords
    ######################################
    pot = distpot.MinPermDistPotential( comcoords1, comcoords2, 0.4, permlist )
    #if False:
        #print some stuff. not necessary
        #Emin = pot.getEnergy(aamin)
        #dist, X11, X22 = findBestPermutationRBMol(coords1, aa2xyz(X2in, aamin), permlist )
        #print "initial energy", Emin, "dist", dist
    saveit = storage.SaveN( 20 )
    takestep = distpot.random_rotation
    bh = basinhopping.BasinHopping( aamin, pot, takestep, storage=saveit.insert)


    Eminglobal = pot.globalEnergyMin() #condition for determining isomer
    print "global Emin", Eminglobal

    ##########################################################################
    # run basin hopping for ninter steps or until the global minimum is found
    # (i.e. determine they are isomers)
    ##########################################################################
    print "using basin hopping to optimize rotations + permutations"
    for i in range(niter):
        bh.run(1)
        Emin = saveit.data[0][0]
        if abs(Emin-Eminglobal) < 1e-6:
            print "isomer found"
            break

    """
    Lower energies generally mean smaller distances, but it's not guaranteed.
    Check a number of the lowest energy structures. To ensure get the correct
    minimum distance structure.
    """
    print "lowest structures found"
    aamin = saveit.data[0][1]
    dmin = 1000.
    for (E, aa) in saveit.data:
        coords2 = coordsApplyRotation(coords2in, aa)
        dist, X11, X22 = findBestPermutationRBMol(coords1, coords2, mysys, permlist )
        print "E %11.5g dist %11.5g" % (E, dist)
        if dist < dmin:
            dmin = dist
            aamin = aa.copy()

    ###################################################################
    #we've optimized the rotation in a permutation independent manner
    #now optimize the permutation
    ###################################################################
    dbefore = getDistaa(coords1, coords2in, mysys)
    coords2 = coordsApplyRotation(coords2in, aamin)
    dafter = getDistaa(coords1, coords2, mysys)
    print "dist before, after applying rotation from basin hopping", dbefore, dafter, mysys.getEnergy(coords2in), mysys.getEnergy(coords2)
    dmin, coords1, coords2min = findBestPermutationRBMol(coords1, coords2, mysys, permlist )
    #dafter = getDistaa(coords1, coords2min, mysys)
    print "dist findBestPerm", dmin, mysys.getEnergy(coords2min)


    ###################################################################
    # permutations are set, do one final mindist improve accuracy
    #of rotation optimization
    ###################################################################
    #dmin, X2min = alignRotation( X1, X2min )

    return dmin, coords1, coords2min

def randomCoords(nmol):
    coords = np.zeros(2*3*nmol, np.float64)
    coords[0:3*nmol] = np.random.uniform(-1,1,[nmol*3]) * 1.3*(3*nmol)**(1./3)
    for i in range(nmol):
        k = 3*nmol + 3*i
        coords[k : k + 3] = rot.random_aa()
    return coords


def test(coords1, coords2, mysys, permlist):
    printlist = []
    printlist.append((coords1.copy(), "after quench"))
    printlist.append((coords2.copy(), "after quench"))
    coords2in = coords2.copy()


    distinit = getDistaa(coords1, coords2, mysys)
    print "distinit", distinit

    (dist, coords1, coords2) = minPermDistRBMol(coords1,coords2, mysys, permlist=permlist)
    distfinal = getDistaa(coords1, coords2, mysys)
    print "dist returned    ", dist
    print "dist from coords ", distfinal
    print "initial energy", mysys.getEnergy(coords2in)
    print "final energy  ", mysys.getEnergy(coords2)

    printlist.append((coords1.copy(), "coords1 after mindist"))
    printlist.append((coords2.copy(), "coords2 after mindist"))
    import printing.print_atoms_xyz as printxyz
    with open("otp.xyz", "w") as fout:
        for coords, line2 in printlist:
            xyz = mysys.getxyz(coords)
            printxyz.printAtomsXYZ(fout, xyz, line2=line2, atom_type = ["N", "O", "O"])
    
    
    return dist, coords1, coords2


def test_LWOTP(nmol = 5):
    from potentials.rigid_bodies.molecule import Molecule, setupLWOTP
    from potentials.rigid_bodies.sandbox import RBSandbox
    from potentials.lj import LJ

    printlist = []
    
    #set up system
    otp = setupLWOTP()
    #set up a list of molecules
    mols = [otp for i in range(nmol)]
    # define the interaction matrix for the system.
    # for LWOTP there is only one atom type, so this is trivial
    lj = LJ()
    interaction_matrix = [[lj]]
    #set up the RBSandbox object
    mysys = RBSandbox(mols, interaction_matrix)
    nsites = mysys.nsites

    permlist = [range(nmol)]

    #define initial coords
    coords1 = randomCoords(nmol)
    printlist.append( (coords1.copy(), "very first"))
    #quench X1
    ret = quench( coords1, mysys.getEnergyGradient)
    coords1 = ret[0]
    
    #define initial coords2
    coords2 = randomCoords(nmol)
    printlist.append( (coords2.copy(), "very first"))
    #quench X2
    ret = quench( coords2, mysys.getEnergyGradient)
    coords2 = ret[0]
    coords2in = coords2.copy()
    
    printlist.append((coords1.copy(), "after quench"))
    printlist.append((coords2.copy(), "after quench"))

    print "******************************"
    print "testing for OTP not isomer"
    print "******************************"
    d, c1new, c2new = test(coords1, coords2, mysys, permlist)
    
    print "******************************"
    print "testing for OTP ISOMER"
    print "******************************"
    coords1 = coords2in.copy()
    coords2 = c2new.copy()
    #try to reverse the permutations and symmetry operations we just applied 
    test(coords1, coords2, mysys, permlist)
        
        

if __name__ == "__main__":
    print "******************************"
    print "testing for OTP"
    print "******************************"
    test_LWOTP(5)
    
    
