import distpot
import numpy as np
import copy
import rotations as rot
from quench import quench
import basinhopping
import storage.savenlowest as storage
from mindistutils import CoMToOrigin, alignRotation, findBestPermutationRBMol, getDistaa


class PermGroup:
    """
    grouplist is a list of allowed permutations.
    
    Each element in grouplist defines a permutation.
        If that permutation is dependent, it is another PermGroup object
        if that permutation is independent, it is a list of permutable atoms
    
    Each element in grouplist is independent of all the other elements
    
    Each element  
    
    for a minimal PermGroup, grouplist is a single list, or a list 
    of lists of fully permutatble atoms."""
    def __init__(self, grouplist):
        #self.atomlist = [] #the atoms involved in this group
        self.grouplist = [] #the permutation groups in this group
        pass

"""
for 3 water molecules.    O   H H 
permlistmol = [ [0,[ [0],[1,2] ]], [3, [ [3],[4,5] ]], [6,[[6], [7,8] ]] ]

PermGroup.grouplist     = [ H2OPermGroup1, H2OPermGroup2, H2OPermGroup3 ]
H2OPermGroup1.grouplist = [ 0, H2PermGroup1 ]  # two types of permutation
 H2PermGroup1.grouplist = [ [1,2] ]
H2OPermGroup1.grouplist = [ 3, H2PermGroup2 ]
 H2PermGroup1.grouplist = [ [4,5] ]
H2OPermGroup1.grouplist = [ 6, H2PermGroup3 ]
 H2PermGroup1.grouplist = [ [7,8] ]


for 3 O2 molecules        O O     

PermGroup.grouplist    = [ O2PermGroup1, O2PermGroup2, O2PermGroup3 ]
O2PermGroup1.grouplist = [ [0,1] ]  #one type of permutation
O2PermGroup2.grouplist = [ [2,3] ]
O2PermGroup3.grouplist = [ [4,5] ]


the phenylalanine example at http://www-wales.ch.cam.ac.uk/OPTIM.doc/node4.html
PermGroup1.grouplist = [ [[1,2], [3,4], [5,6], [7,8]] ] #1 permutation of 4 dependent molecules
PermGroup1.grouplist = [ [[1,2], [3,4], [5,6], [7,8]] ] #1 permutation of 4 dependent molecules



"""

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
            aamin = aa

    ###################################################################
    #we've optimized the rotation in a permutation independent manner
    #now optimize the permutation
    ###################################################################
    dbefore = getDistaa(coords1, coords2in, mysys)
    coords2 = coordsApplyRotation(coords2in, aamin)
    dafter = getDistaa(coords1, coords2, mysys)
    print "dist before, after applying rotation from basin hopping", dbefore, dafter
    dmin, coords1, coords2min = findBestPermutationRBMol(coords1, coords2, mysys, permlist )
    #dafter = getDistaa(coords1, coords2min, mysys)
    print "dist findBestPerm", dmin


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
    
    printlist.append((coords1.copy(), "after quench"))
    printlist.append((coords2.copy(), "after quench"))


    #X1 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 0., 1.,] )
    #X2 = np.array( [ 0., 0., 0., 1., 0., 0., 0., 1., 0.,] )
    coords1i = copy.copy(coords1)
    coords2i = copy.copy(coords2)

    distinit = np.linalg.norm(coords1-coords2)
    print "distinit", distinit

    (dist, coords1, coords2) = minPermDistRBMol(coords1,coords2, mysys, permlist=permlist)
    distfinal = getDistaa(coords1, coords2, mysys)
    print "dist returned    ", dist
    print "dist from coords ", distfinal


    printlist.append((coords1.copy(), "coords1 after mindist"))
    printlist.append((coords2.copy(), "coords2 after mindist"))
    import printing.print_atoms_xyz as printxyz
    with open("otp.xyz", "w") as fout:
        for coords, line2 in printlist:
            xyz = mysys.getxyz(coords)
            printxyz.printAtomsXYZ(fout, xyz, line2=line2, atom_type = ["N", "O", "O"])
        
        

if __name__ == "__main__":
    print "******************************"
    print "testing for OTP"
    print "******************************"
    test_LWOTP(5)
    
    
