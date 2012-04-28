import numpy as np
import molecule
import rotations as rot


class RBSandbox:
    """
    a system of molecules
    """
    def __init__(self, molecule_list):
        """
        molecule list is a Molecule objects that define the system
        """
        self.molecules = molecule_list
        self.nmol = len(self.molecules)
        self.nsites = 0
        for mol in self.molecules:
            self.nsites += mol.nsites
        
    def getxyz(self, comcoords, aacoords):
        print len(comcoords), 3*self.nmol
        print len(aacoords), 3*self.nmol
        assert(len(comcoords) == 3*self.nmol)
        assert(len(aacoords) == 3*self.nmol)
        xyz = np.zeros(self.nsites*3)
        k = 0
        for imol, mol in enumerate(self.molecules):
            molxyz = mol.getxyz( aacoords[imol*3:imol*3+3] )
            print "before ", molxyz
            for j in range(mol.nsites):
                molxyz[j*3:j*3+3] += comcoords[imol*3:imol*3+3]
            print "after  ", molxyz

            xyz[k:k+mol.nsites*3] = molxyz
            print "xyz", xyz
            k += mol.nsites*3
        return xyz
            
            


def test_molecule():
    from numpy import sin, cos, pi
    import copy
    otp = molecule.setupLWOTP()
        
    xyz = otp.getxyz()
    from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    import sys
    #with open("out.xyz", "w") as fout:
    printxyz(sys.stdout, xyz)
    
    nmol = 4
    mols = [otp for i in range(nmol)]
    mysys = RBSandbox(mols)
    comcoords = np.random.uniform(-1,1,[nmol*3]) * 2*(nmol)**(1./3)
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
    
    xyz = mysys.getxyz( coords[0:3*nmol], coords[3*nmol:] )
    
    printlist = []
    printlist.append( (xyz.copy(), "initial"))
    
    #from printing.print_atoms_xyz import printAtomsXYZ as printxyz
    with open("otp.xyz", "w") as fout:
        for xyz, line2 in printlist:
            printxyz( fout, xyz, line2=line2)
    
    
    
if __name__ == "__main__":
    test_molecule()