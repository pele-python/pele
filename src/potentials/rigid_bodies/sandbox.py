import numpy as np
import molecule
import rotations as rot
import itertools
from potentials.potential import potential
import copy


class RBSandbox(potential):
    """
    Defines a system of rigid body molecules
    """
    def __init__(self, molecule_list, interaction_matrix):
        """
        molecule_list:  a list of Molecule objects that define the system

        interaction_matrix:  a matrix that defines the interactions between site types.
            inter = interaction_matrix[site1.type][site2.type] is the interaction between site1 and site2
            inter.getEnergy(coords) returns the energy of the interaction
        """
        self.molecule_list = [copy.deepcopy(mol) for mol in molecule_list]
        self.interaction_matrix = interaction_matrix
        self.nmol = len(self.molecule_list)

        self.nsites = 0
        self.typelist = []
        self.nsites_cum = np.zeros(self.nmol, np.int)
        for i, mol in enumerate(self.molecule_list):
            self.nsites_cum[i] = self.nsites
            self.nsites += mol.nsites
            for site in mol.sitelist:
                self.typelist.append( site.type )

        self.oldcoords = np.zeros(3*2*self.nmol)

    def getxyz(self, coords):
        """convert center of mass + angle-axis coords into xyz coordinates"""
        self.update(coords)
        xyz = np.zeros(self.nsites*3)
        isite = 0
        for mol in self.molecule_list:
            for site in mol.sitelist:
                xyz[isite*3 : isite*3 + 3] = site.abs_position
                isite += 1
        return xyz

    def molmolEnergyOld(self, imol1, imol2, xyz):
        E = 0.
        mol1 = self.molecule_list[imol1]
        mol2 = self.molecule_list[imol2]
        #print "mol ", imol1, imol2
        k1 = self.nsites_cum[imol1]*3
        k2 = self.nsites_cum[imol2]*3
        coords2 = np.zeros(6, np.float64)
        for i1, site1 in enumerate(mol1.sitelist):
            coords2[0:3] = xyz[k1 + i1*3 : k1 + i1*3+3]
            for i2, site2 in enumerate(mol2.sitelist):
                interaction = self.interaction_matrix[site1.type][site2.type]
                coords2[3:] = xyz[k2 + i2*3 : k2 + i2*3+3]
                dE = interaction.getEnergy(coords2)
                E += dE
                #dr = np.linalg.norm(coords2[0:3] - coords2[3:])
                #print "    site", k1 + i1*3, k2 + i2*3, dE, dr

        return E

    def molmolEnergy(self, mol1, mol2):
        """
        return the energy of interaction between mol1 and mol2.

        don't update the molecules
        """
        E = 0.
        coords2 = np.zeros(6, np.float64)
        for site1 in mol1.sitelist:
            coords2[0:3] = site1.abs_position
            for site2 in mol2.sitelist:
                coords2[3:6] = site2.abs_position
                #print coords2
                interaction = self.interaction_matrix[site1.type][site2.type]
                dE = interaction.getEnergy(coords2)
                E += dE
        return E


    def getEnergy(self, coords):
        self.update(coords)
        E = 0.
        for imol1 in range(self.nmol):
            mol1 = self.molecule_list[imol1]
            for imol2 in range(0,imol1):
                mol2 = self.molecule_list[imol2]
                E += self.molmolEnergy(mol1, mol2)
        return E



    def molmolEnergyGradient(self, mol1, mol2):
        """
        Update the energy and gradients on mol1 and mol2
        """
        E = 0.
        coords2 = np.zeros(6, np.float64)
        for i1, site1 in enumerate(mol1.sitelist):
            coords2[0:3] = site1.abs_position
            for i2, site2 in enumerate(mol2.sitelist):
                coords2[3:6] = site2.abs_position
                #print coords2
                interaction = self.interaction_matrix[site1.type][site2.type]
                dE, V2 = interaction.getEnergyGradient(coords2)
                E += dE
                mol1.E += dE
                mol2.E += dE
                mol1.comgrad += V2[0:3]
                mol2.comgrad += V2[3:6]
                #do angle axis part
                #calculate drdp
                mol1.aagrad += np.dot( site1.drdp, V2[0:3])
                mol2.aagrad += np.dot( site2.drdp, V2[3:6])

    def coords_compare(self, coords):
        """ return true if coords is the same as oldcoords"""
        maxdif = np.max( np.abs(coords - self.oldcoords))
        return maxdif < 1e-16
            #print "maxdif", maxdif, coords[0:3]
            #return True
        #return False

    def update(self, coords):
        """
        update the com and angle-axis coords and dependents on all molecules

        only do it if coords is different.
        """
        #using coords_compare makes quench almost always fail for some reason.
        #if self.coords_compare( coords):
            #return
        self.oldcoords[:] = coords
        nmol= self.nmol
        for imol, mol in enumerate(self.molecule_list):
            com = coords[         imol*3 :          imol*3 + 3]
            aa  = coords[3*nmol + imol*3 : 3*nmol + imol*3 + 3]
            mol.update( com, aa )
        #for imol, mol in enumerate(self.molecule_list):
            #print "abs position", imol, mol.com


    def getEnergyGradient(self, coords):
        """
        return energy and gradient

        update energy and gradient on all the molecules
        """
        nmol = self.nmol
        self.update(coords)
        Etot = 0.
        grad = np.zeros([2*self.nmol*3])
        for imol1 in range(self.nmol):
            mol1 = self.molecule_list[imol1]
            for imol2 in range(imol1):
                #print "mol mol", imol1, imol2
                mol2 = self.molecule_list[imol2]
                self.molmolEnergyGradient(mol1, mol2)
        for imol, mol in enumerate(self.molecule_list):
            Etot += mol.E
            grad[         3*imol :          3*imol + 3] = mol.comgrad
            #print "aagrad", mol.aagrad
            grad[3*nmol + 3*imol : 3*nmol + 3*imol + 3] = mol.aagrad
        Etot /= 2.
        return Etot, grad


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



def test_sandbox(nmol = 6):
    from numpy import sin, cos, pi
    import copy
    from potentials.lj import LJ

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
    import quench
    coords, E, rms, funcalls = quench.quench(coords, mysys.getEnergyGradient, iprint=-1)
    print "postquench E", E, "rms", rms, "funtion calls", funcalls
    xyz = mysys.getxyz(coords )
    printlist.append( (xyz.copy(), "post quench"))


    #print the saved coords
    fname = "otp.xyz"
    print "saving xyz coords to", fname
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    with open(fname, "w") as fout:
        for xyz, line2 in printlist:
            printxyz( fout, xyz, line2=line2, atom_type=["N", "O", "O"])
            
    
    test_symmetries(coords, mysys)



if __name__ == "__main__":
    test_sandbox()

