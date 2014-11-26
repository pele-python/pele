from playground.molecule.molecule import Molecule, Atom
from playground.molecule.parser import PymolParser
from playground.molecule.molecularsystem import MolecularSystem
import pele.utils.elements as elem 

import networkx as nx
import unittest
import os

class TestMolecule(unittest.TestCase):
    """
    a base class for molecule unit tests
    """

    def setUp(self):
        ''' Function Creates a reference molecule system with known data'''
        # Generate comparison coords list
        self.coords = [4.954, -0.924, -5.684,
                       5.427, 0.193, -4.880,
                       5.873, -0.251, -3.503,
                       5.756, -1.417, -3.054]

        self.atom_names = ['N', 'S', 'CA', 'O']
        self.atom_symbols = [elem.alternate_names[name] for name in self.atom_names] 

        # Generate the Atom objects - zero based index ordered as in file.
        N = Atom(0, 'N')
        S = Atom(1, 'S')
        CA = Atom(2, 'CA')
        O = Atom(3, 'O')

        # Create the graph
        self.mol_graph = nx.Graph()

        # add the nodes in the order of appearance in the file.
        self.mol_graph.add_node(N)
        self.mol_graph.add_node(S)
        self.mol_graph.add_node(CA)
        self.mol_graph.add_node(O)

        # add the edges. don't care so much about the order here.
        self.mol_graph.add_edge(N, S)
        self.mol_graph.add_edge(S, CA)
        self.mol_graph.add_edge(CA, O)

        # create a molecule system from the test data
        reference_molecule = Molecule(0, self.coords, self.mol_graph)

        # create a molecular system from the manual data
        self.reference_molecular_system = MolecularSystem([reference_molecule])

    def test_parse_xyz(self):
        # Generate test data
        xyz_data = []
        num_atoms = len(self.atom_names)
        xyz_data.append("    " + str(num_atoms) + "\n")
        xyz_data.append("    Energy Comment\n")
        
        coord_index = 0
        for atom_symbol in self.atom_symbols:
            l='{:1} {: >8.3f}{: >8.3f}{: >8.3f}\n'.format(atom_symbol, 
                                                          self.coords[3 * coord_index],
                                                          self.coords[3 * coord_index + 1],
                                                          self.coords[3 * coord_index + 2])
            xyz_data.append(l)
            coord_index+=1
            
        # Save test data as a file
        try:
            with open('test.xyz', "w") as f_out:
                f_out.writelines(xyz_data)
        except IOError:
            raise Exception("Unable to write xyz file for output")

        self.check_Molecule_System('test.xyz', 'pymol')

        # delete the test file
        os.remove('test.xyz')

    def test_molecule_pdb(self):
        # Generate test data
        pdb_data = []

        coord_index = 0
        for atom_name, atom_symbol in zip(self.atom_names, self.atom_symbols):
            l='ATOM {: >06d} {: <4} GLY A   1    {: >8.3f}{: >8.3f}{: >8.3f}  1.00  5.04           {: >2}\n'.format(coord_index+1,
                                                                                                                       atom_name,
                                                                                                                       self.coords[3 * coord_index],
                                                                                                                       self.coords[3 * coord_index + 1],
                                                                                                                       self.coords[3 * coord_index + 2],
                                                                                                                       atom_symbol)
            pdb_data.append(l)
            coord_index+=1

        # Save test data as a file
        try:
            with open('test.pdb', "w") as f_out:
                f_out.writelines(pdb_data)
        except IOError:
            raise Exception("Unable to write pdb file for output")

        self.check_Molecule_System('test.pdb', 'pymol')

        # delete the test file
        os.remove('test.pdb')

    def check_Molecule_System(self, filename, parser):

        # #check for the specified parser
        if parser == 'pymol':
            # Create the Molecule class for testing
            parser = PymolParser(filename, '-qc')
            test_molecular_system = parser.get_molecular_system()
        
        # confirm that test molecule representation is identical to
        # the molecule created locally.

        # ensure that both networks are isomorphic - same number of
        # nodes and same edge connection pattern.
        l_isomorphic = nx.is_isomorphic(self.reference_molecular_system.molecules[0].topology, test_molecular_system.molecules[0].topology)
        self.assertTrue(l_isomorphic)

        # if the network is isomorphic check that the coords and labels of the nodes correspond in
        # both the test and reference molecule.
        if l_isomorphic:
            
            # construct dictionaries of the atomic coordinates keyed by the symbols 
            coords_labels_ref = {}
            coords_labels_test = {}
            
            # populate the dictionaries by looking at each node and extracting the relevant coords from the coords array.
            for atom in test_molecular_system.molecules[0].topology.nodes():
                coords_labels_test[atom.symbol] = (test_molecular_system.molecules[0].coords[ 3 * atom.id],
                                                   test_molecular_system.molecules[0].coords[ 3 * atom.id + 1],
                                                   test_molecular_system.molecules[0].coords[ 3 * atom.id + 2])

            for atom in self.reference_molecular_system.molecules[0].topology.nodes():
                coords_labels_ref[atom.symbol] = (self.reference_molecular_system.molecules[0].coords[ 3 * atom.id], 
                                                  self.reference_molecular_system.molecules[0].coords[ 3 * atom.id + 1], 
                                                  self.reference_molecular_system.molecules[0].coords[ 3 * atom.id + 2])
            
            # loop through the atom names and check they index the same coords.
            for label in self.atom_symbols:
                ref_coords = coords_labels_ref[label]
                test_coords = coords_labels_test[label]
                for ref_val, test_val in zip(ref_coords, test_coords):
                    self.assertAlmostEqual(ref_val,test_val, delta = 0.00001)
                
            ## Old test. doesn't work since we may renumber the atoms completely in the parser class but nifty code.  
            # ensure the edges in the test molecule connect the same atoms
            # as the edges in the constructed molecule
            # bond_labels = {tuple(sorted([bond[0].id, bond[1].id])) for bond in self.mol_graph.edges_iter()}
            # test_bond_labels = {tuple(sorted([test_bond[0].id, test_bond[1].id])) for test_bond in test_molecular_system.molecules[0].topology.edges_iter()}

            # check that the sets are equivalent (a permutation of each other)
            # self.assertTrue(test_bond_labels == bond_labels)

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
