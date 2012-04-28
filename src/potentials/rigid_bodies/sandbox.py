import numpy as np
import molecule
import rotations as rot
import itertools
from potentials.potential import potential


class RBSandbox(potential):
    """
    a system of molecule_list
    """
    def __init__(self, molecule_list, interaction_matrix):
        """
        molecule_list:  a list of Molecule objects that define the system
        
        interaction_matrix:  a matrix that defines the interactions between site types.
            inter = interaction_matrix[site1.type][site2.type] is the interaction between site1 and site2
            inter.getEnergy(coords) returns the energy of the interaction
        """
        self.molecule_list = molecule_list
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
        
    def getxyz(self, coords):
        nmol = self.nmol
        xyz = np.zeros(self.nsites*3)
        for imol, mol in enumerate(self.molecule_list):
            molxyz = mol.getxyz( coords[3*nmol + imol*3 : 3*nmol + imol*3+3] )
            for j in range(mol.nsites):
                molxyz[j*3:j*3+3] += coords[imol*3 : imol*3+3]
            
            k = self.nsites_cum[imol]*3
            xyz[k:k+mol.nsites*3] = molxyz
        return xyz
    
    def molmolEnergy(self, imol1, imol2, xyz):
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
                dr = np.linalg.norm(coords2[0:3] - coords2[3:])
                #print "    site", k1 + i1*3, k2 + i2*3, dE, dr

        return E

    
    def getEnergy(self, coords):
        xyz = self.getxyz(coords)
        E = 0.
        for imol1 in range(self.nmol):
            for imol2 in range(0,imol1):
                E += self.molmolEnergy(imol1, imol2, xyz)            
        return E


def test_sandbox():
    from numpy import sin, cos, pi
    import copy
    from potentials.lj import LJ
    lj = LJ()
    interaction_matrix = [[lj]]

    
    otp = molecule.setupLWOTP()
        
    xyz = otp.getxyz()
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    import sys
    #with open("out.xyz", "w") as fout:
    printxyz(sys.stdout, xyz)
    
    nmol = 4
    mols = [otp for i in range(nmol)]
    mysys = RBSandbox(mols, interaction_matrix)
    nsites = mysys.nsites
    comcoords = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
    aacoords = np.array( [copy.copy(rot.random_aa()) for i in range(nmol)] )
    print comcoords
    print aacoords
    aacoords = aacoords.reshape(3*nmol)
    print aacoords
    coords = np.zeros(2*3*nmol, np.float64)
    print "lencoords, len aacoords", len (coords), len(aacoords), len(comcoords)

    coords[0:3*nmol] = comcoords[:]
    coords[3*nmol:2*3*nmol] = aacoords[:]
    print "lencoords, len aacoords", len (coords), len(aacoords), len(comcoords)
    
    xyz = mysys.getxyz( coords )
    
    printlist = []
    printlist.append( (xyz.copy(), "initial"))
    
            
    E = mysys.getEnergy(coords)
    print "energy", E
    
    numericV = mysys.NumericalDerivative(coords, 1e-4)
    print "numeric V", numericV
    
    import quench
    coords, E, rms, funcalls = quench.quench(coords, mysys.getEnergyGradient, iprint=1)
    print "postquench", E, rms, funcalls
    print coords
    xyz = mysys.getxyz(coords )
    printlist.append( (xyz.copy(), "post quench"))

    
    #from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    with open("otp.xyz", "w") as fout:
        for xyz, line2 in printlist:
            printxyz( fout, xyz, line2=line2)

    
    
if __name__ == "__main__":
    test_sandbox()
    