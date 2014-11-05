from molecule import Molecule, Atom, Bond
from utilities import write_text_file,read_text_file


import sys
import os
import urllib


class TestMolecule(unittest.TestCase):
    """
    a base class for molecule unit tests
    """
    def runtest(self, X1, X2, mindist):
        X1i = np.copy(X1)
        X2i = np.copy(X2)
        
        (distreturned, X1, X2) = mindist(X1, X2)

        distinit = np.linalg.norm(self.X1 - X2i)
        distfinal = np.linalg.norm(X1 - X2)
        self.assertAlmostEqual( distfinal, distreturned, 5, "returned distance is wrong: %g != %g, difference %g" % (distfinal, distreturned, distfinal-distreturned) )
        self.assertLessEqual(distfinal, distinit)
        
        #test if the energies have changed
        Ei = self.pot.getEnergy(X1i)        
        Ef = self.pot.getEnergy(X1)
        self.assertAlmostEqual( Ei, Ef, 10, "Energy of X1 changed: %g - %g = %g" % (Ei, Ef, Ei - Ef) )
        Ei = self.pot.getEnergy(X2i)        
        Ef = self.pot.getEnergy(X2)
        self.assertAlmostEqual( Ei, Ef, 10, "Energy of X2 changed: %g - %g = %g" % (Ei, Ef, Ei - Ef) )

        return distreturned, X1, X2

def test_read_write():
    lines=["The lazy brown dog jumped over the sly red fox\n","The lazy brown dog jumped over the sly red fox\n"]
    write_text_file(lines,'testfile.txt')
    newlines=read_text_file('testfile.txt')
    os.remove('testfile.txt')
    
    return not False in [ False for newline,line in zip(newlines,lines) if newline!=line]

def test_pymol():
        #Set up Pymol
        pymol=pymol_wrapper.init('-qc')
            
        #load the data into pymol    
        self.load_data_pymol(pymol,l_input_filename)

        #Extract the pymol object defining the molecule
        l_pymol_data=self.extract_data_pymol(pymol)


def test_extract_atoms():
    
def test_extract_data():
    
def test_extract_topology():
    

def test_molecule:


def test_atom():
    atom=Atom(1,'H')
    result=False
    if atom.atom_id==1and atom.atom_name=='H':
        result=True
    return result
    
def test_bond():
    bond=Bond(1,1,2,1)
    result=False
    if bond.bond_id==1 and bond.atom1_id=='H' and bond.atom2_id==2 and bond.order==1:
        result=True
    return result

def fetch_pdb(id):
  url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
  return urllib.urlopen(url).read()


if __name__=='__main__':

    testResults=[]
    testResults.append(test_atom())
    testResults.append(test_bond())

    #create a "molecule"
    test_data_coords=[0.0 0.0 0.0 
                      0.0 1.0 0.0 
                      0.0 1.0 -1.0]
    
    atom1=Atom(1,'H')
    atom1=Atom(2,'O')
    atom1=Atom(3,'H')
    
    

    
    
    testResults.append(test_read_write())
    
    
    
    testStructure = write_text_file(fetch_pdb('2RNM'),'2RNM.pdb')
    
    
    HETS=Molecule('2RNM.pdb')
    
    HETS.coords
    HETS.extract_atoms(l_pymol_object_data)
    HETS.extract_data_pymol(pymol)
    HETS.extract_topology(l_pymol_object_data, l_pymol_index_map)
    HETS.load_data_pymol(pymol, l_input_filename)
    HETS.topology
    HETS.__init__(l_input_filename)
    
    
    
    
    
