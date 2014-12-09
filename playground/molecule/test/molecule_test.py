from playground.molecule.molecule import Molecule, Atom
from playground.molecule.molecularsystem import MolecularSystem
import pele.utils.elements as elem

import networkx as nx
import numpy as np
import unittest
import os


class TestMolecule(unittest.TestCase):
    """
    a base class for molecular system unit tests
    """

    def setUp(self):
        """ Function creates a reference molecular system with known data"""
        
        self.generate_ref_system('basic');
        
    def generate_ref_system(self, ref_system_type):
        
        # generate a reference system that depends on the keyword selected
        if ref_system_type == 'basic':
            self.gen_basic_ref_system()
        else:
            # default to basic refer system
            self.gen_basic_ref_system()
        
    def gen_basic_ref_system(self):
        
        # Generate comparison coords list
        self.coords = [4.954, -0.924, -5.684,
                       5.427, 0.193, -4.880,
                       5.873, -0.251, -3.503,
                       5.756, -1.417, -3.054]

        self.atom_names = ['N', 'S', 'CA', 'O']
        self.atom_symbols = [elem.alternate_names[name] for name in self.atom_names]

        self.residue = 'GLY'
        self.resid = 1
        self.chain = 'A'

        # Generate the Atom objects - zero based index ordered as in file.
        N = Atom(0, 'N')
        S = Atom(1, 'S')
        CA = Atom(2, 'CA')
        O = Atom(3, 'O')

        # Create the graph
        self.mol_graph = nx.Graph()

        # add the nodes along with an attribute - the symbol - used to test equivalence of nodes between molcules.
        self.mol_graph.add_node(N, atom=N)
        self.mol_graph.add_node(S, atom=S)
        self.mol_graph.add_node(CA, atom=CA)
        self.mol_graph.add_node(O, atom=O)

        # add the edges. don't care so much about the order here.
        self.mol_graph.add_edge(N, S, hweight = np.sqrt(N.mass * S.mass))
        self.mol_graph.add_edge(S, CA, hweight = np.sqrt(S.mass * CA.mass))
        self.mol_graph.add_edge(CA, O, hweight = np.sqrt(CA.mass * O.mass))

        # create a molecule system from the test data
        reference_molecule = Molecule(0, self.coords, self.mol_graph)

        # create a molecular system from the data
        self.reference_molecular_system = MolecularSystem()
        self.reference_molecular_system.add_molecule(reference_molecule)

    def test_parse_xyz(self):
        # Generate test data - dumps the reference system to an xyz file.
        xyz_data = []
        num_atoms = len(self.atom_names)
        xyz_data.append("    " + str(num_atoms) + "\n")
        xyz_data.append("    Energy Comment\n")

        coord_index = 0
        for atom_symbol in self.atom_symbols:
            l = '{:1} {: >8.3f}{: >8.3f}{: >8.3f}\n'.format(atom_symbol,
                                                            self.coords[3 * coord_index],
                                                            self.coords[3 * coord_index + 1],
                                                            self.coords[3 * coord_index + 2])
            xyz_data.append(l)
            coord_index += 1

        # Save test data as a file
        try:
            with open('test.xyz', "w") as f_out:
                f_out.writelines(xyz_data)
        except IOError:
            raise Exception("Unable to write xyz file for output")

        # calls the testing function
        self.check_Molecule_System('test.xyz', 'pymol', args='-qc')

        # delete the test file
        os.remove('test.xyz')

    def test_molecule_pdb(self):
        # Generate test data -  dumps the reference data to a dummy PDB
        pdb_data = []

        coord_index = 0
        for atom_name, atom_symbol in zip(self.atom_names, self.atom_symbols):
            l = 'ATOM {: >06d} {: <4} {:3} {:1}{: >4d}    {: >8.3f}{: >8.3f}{: >8.3f}  1.00  5.04           {: >2}\n'.format(
                coord_index + 1,
                atom_name,
                self.residue,
                self.chain,
                int(self.resid),
                self.coords[3 * coord_index],
                self.coords[3 * coord_index + 1],
                self.coords[3 * coord_index + 2],
                atom_symbol)
            pdb_data.append(l)
            coord_index += 1

        # Save test data as a file
        try:
            with open('test.pdb', "w") as f_out:
                f_out.writelines(pdb_data)
        except IOError:
            raise Exception("Unable to write pdb file for output")

        # calls the checking function to assert that the test molecular system
        # and reference molecular system are the same.

        self.check_Molecule_System('test.pdb', parser='pymol', args='-qc')

        # delete the test file
        os.remove('test.pdb')

    def check_Molecule_System(self, l_filename, parser='pymol', args='-qc'):

        # Create the molecular system
        self.test_molecular_system = MolecularSystem()
        
        # use the molecular system to load in the filename and create the molecule.
        self.test_molecular_system.load_file(l_filename, parser, args)
        
        # Check that the molecular system loaded into the test_molecular system
        # is equal to the molecular system in the reference molecular system 
        self.assertEqual(self.test_molecular_system, self.reference_molecular_system)
        

    def test_atom(self):
        atom = Atom(1, '1H')
        self.assertEqual(1, atom.id)
        self.assertEqual('1H', atom.alt_symbol)
        self.assertEqual('hydrogen', atom.name)
        self.assertEqual('H', atom.symbol)
        self.assertEqual([1.000, 1.000, 1.000], atom.color)
        self.assertEqual(1.00794000, atom.mass)
        self.assertEqual(1.2000, atom.radius)

        atom2 = Atom(2, 'Na')
        self.assertEqual(2, atom2.id)
        self.assertEqual('Na', atom2.alt_symbol)
        self.assertEqual('sodium', atom2.name)
        self.assertEqual('Na', atom2.symbol)
        self.assertEqual([0.671, 0.361, 0.949], atom2.color)
        self.assertEqual(22.98976928, atom2.mass)
        self.assertEqual(2.2700, atom2.radius)


if __name__ == '__main__':
    unittest.main()
