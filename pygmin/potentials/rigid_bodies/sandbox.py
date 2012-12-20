import numpy as np
import molecule
import pygmin.utils.rotations as rot
#from potentials.rigid_body_potential import RigidBodyPotential
from rigid_body_system import RigidBodySystem

class RBSandbox(RigidBodySystem):
    def __init__(self, molecule_list, interaction_matrix):
        """
        a wrapper class to combine RigidBodySystem with RBSandboxPotential

        molecule_list:  a list of Molecule objects that define the system

        interaction_matrix:  a matrix that defines the interactions between site types.
            inter = interaction_matrix[site1.type][site2.type] is the interaction between site1 and site2
            inter.getEnergy(coords) returns the energy of the interaction
        """
        RigidBodySystem.__init__(self, molecule_list)
        from sandbox_potential import RBSandboxPotential
        self.potential = RBSandboxPotential(self, interaction_matrix)
        #self.setPotential( self.potential )



def test_symmetries(coords, mysys):
    nmol = mysys.nmol
    #the molecules now have their correct permutation.  
    #For each molecule, apply the symmetry operation which minimized the distance
    for i, mol in enumerate(mysys.molecule_list):
        comold = coords[3*i : 3*i + 3]
        kaa = 3*nmol + 3*i
        aaold = coords[kaa : kaa + 3]
        for xyz, aanew in mol.getSymmetries(comold, aaold):
            coords[kaa : kaa + 3] = aanew
            print "E symmetry", mysys.getEnergy(coords) 


def random_coords(nmol, nsites = None):
    if nsites == None: nsites = nmol
    coords = np.zeros(2*3*nmol)
    coords[0:nmol*3] = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    for i in range(nmol):
        k = nmol*3 + 3*i
        coords[k : k + 3] = rot.random_aa()
    return coords


def test_sandbox(nmol = 6):
    import copy
    from pygmin.potentials.lj import LJ

    #define the molecule types.
    #here use only one type, LWOTP
    otp = molecule.setupLWOTP()


    # define the interaction matrix for the system.
    # for LWOTP there is only one atom type, so this is trivial
    lj = LJ()
    interaction_matrix = [[lj]]


    #set up a list of molecules
    mols = [otp for i in range(nmol)]

    #set up the RBSandbox object
    mysys = RBSandbox(mols, interaction_matrix)
    nsites = mysys.nsites
    

    #get an initial set of coordinates
    comcoords = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    aacoords = np.array( [copy.copy(rot.random_aa()) for i in range(nmol)] )
    aacoords = aacoords.reshape(3*nmol)
    coords = np.zeros(2*3*nmol, np.float64)
    coords[0:3*nmol] = comcoords[:]
    coords[3*nmol:2*3*nmol] = aacoords[:]
    print "lencoords, len aacoords", len (coords), len(aacoords), len(comcoords)

    #print "xyz coords", mysys.transformToXYZ(coords)

    #save the initial set of coords
    printlist = []
    xyz = mysys.getxyz( coords )
    printlist.append( (xyz.copy(), "initial"))

    #calculate the initial energy
    Einit = mysys.getEnergy(coords)
    print "initial energy", Einit

    #test the gradient
    numericV = mysys.NumericalDerivative(coords, 1e-10)
    numericV = numericV.copy()
    print "numeric V", numericV

    E, V = mysys.getEnergyGradient(coords)
    print "energy from gradient", E, "difference", E - Einit
    #print "analytic gradient", V
    maxgrad_relative = np.max(np.abs(V-numericV)/np.abs(V))
    maxgraddiff = np.max(np.abs(V - numericV))
    print "max error in gradient", maxgraddiff, "max relative", maxgrad_relative

    #do a quench to make sure everything is working well
    from pygmin.optimize import as lbfgs_scipy
    coords, E, rms, funcalls = lbfgs_scipy(coords, mysys.getEnergyGradient, iprint=-1)
    print "postquench E", E, "rms", rms, "funtion calls", funcalls
    xyz = mysys.getxyz(coords )
    printlist.append( (xyz.copy(), "post quench"))


    #print the saved coords
    fname = "otp.xyz"
    print "saving xyz coords to", fname
    from pygmin.printing.print_atoms_xyz import printAtomsXYZ as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz( fout, xyz, line2=line2, atom_type=["N", "O", "O"])
            
    
    test_symmetries(coords, mysys)

def test_sandbox_dumbbell():
    print "***********************************"
    print "testing dumbbell molecule"
    print "***********************************"

    dbel, imatrix = molecule.dumbbell()
    nmol = 5
    print imatrix[0][0] == imatrix[0][1]
    print imatrix[1][0] == imatrix[0][1]
    
    #set up a list of molecules
    mols = [dbel for i in range(nmol)]

    #set up the RBSandbox object
    mysys = RBSandbox(mols, imatrix)
    nsites = mysys.nsites

    coords = random_coords(nmol, nsites)
    
    #calculate the initial energy
    Einit = mysys.getEnergy(coords)
    print "initial energy", Einit

    #test the gradient
    numericV = mysys.NumericalDerivative(coords, 1e-10)
    numericV = numericV.copy()
    print "numeric V", numericV

    E, V = mysys.getEnergyGradient(coords)
    print "energy from gradient", E, "difference", E - Einit
    #print "analytic gradient", V
    maxgrad_relative = np.max(np.abs(V-numericV)/np.abs(V))
    maxgraddiff = np.max(np.abs(V - numericV))
    print "max error in gradient", maxgraddiff, "max relative", maxgrad_relative

    #do a quench to make sure everything is working well
    from pygmin.optimize import lbfgs_scipy
    coords, E, rms, funcalls =  lbfgs_scipy(coords, mysys.getEnergyGradient, iprint=-1)
    print "postquench E", E, "rms", rms, "funtion calls", funcalls

if __name__ == "__main__":
    test_sandbox()
    test_sandbox_dumbbell()

