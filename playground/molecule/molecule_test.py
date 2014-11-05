from molecule import Molecule, Atom, Bond

import networkx as nx
import unittest
import sys
import os
import urllib


class TestMolecule(unittest.TestCase):
    """
    a base class for molecule unit tests
    """

    def test_molecule_xyz(self):
        # Generate test data
        pdb_data = []
        pdb_data.append("    4\n")
        pdb_data.append("    Energy Comment\n")
        pdb_data.append("N   4.954  -0.924  -5.684\n")
        pdb_data.append("C   5.427   0.193  -4.880\n")
        pdb_data.append("C   5.873  -0.251  -3.503\n")
        pdb_data.append("O   5.756  -1.417  -3.054\n")

        # Save test data as a file
        try:
            with open('test.xyz', "w") as f_out:
                f_out.writelines(pdb_data)
        except IOError:
            raise Exception("Unable to write xyz file for output")

        self.check_Molecule_Init('test.xyz', '-qc')

        # delete the test file
        os.remove('test.xyz')

    def test_molecule_pdb(self):
        # Generate test data
        pdb_data = []
        pdb_data.append("ATOM      1  N   GLY A   1       4.954  -0.924  -5.684  1.00  4.90           N\n")
        pdb_data.append("ATOM      2  CA  GLY A   1       5.427   0.193  -4.880  1.00  4.21           C\n")
        pdb_data.append("ATOM      3  C   GLY A   1       5.873  -0.251  -3.503  1.00  4.62           C\n")
        pdb_data.append("ATOM      4  O   GLY A   1       5.756  -1.417  -3.054  1.00  5.04           O\n")

        # Save test data as a file
        try:
            with open('test.pdb', "w") as f_out:
                f_out.writelines(pdb_data)
        except IOError:
            raise Exception("Unable to write pdb file for output")

        self.check_Molecule_Init('test.pdb', '-qc')

        # delete the test file
        os.remove('test.pdb')

    def check_Molecule_Init(self, filename, args):

        # Create the Molecule class for testing
        test_molecule = Molecule(filename, args)

        # Generate comparison coords list
        coords = [4.954, -0.924, -5.684,
                  5.427, 0.193, -4.880,
                  5.873, -0.251, -3.503,
                  5.756, -1.417, -3.054]

        # Generate the Atom objects - zero based index ordered as in file.
        N = Atom(0, 'N')
        CA = Atom(1, 'CA')
        C = Atom(2, 'C')
        O = Atom(3, 'O')

        # Create the bond objects using atoms ids.
        bond_N_CA = Bond(1, 0, 1, 1)
        bond_CA_C = Bond(2, 1, 2, 1)
        bond_C_O = Bond(3, 2, 3, 2)

        # Create the graph
        mol_graph = nx.Graph()

        # add the nodes in the order of appearance in the file.
        mol_graph.add_node(N)
        mol_graph.add_node(CA)
        mol_graph.add_node(C)
        mol_graph.add_node(O)

        # add the edges. don't care so much about the order here.
        mol_graph.add_edge(N, CA, object=bond_N_CA)
        mol_graph.add_edge(CA, C, object=bond_CA_C)
        mol_graph.add_edge(C, O, object=bond_C_O)

        # confirm that test molecule representation is identical to
        # the molecule created locally.

        # Check the coords are the same - same as the order in file.
        for coord, test_coord in zip(coords, test_molecule.coords):
            self.assertAlmostEqual(coord, test_coord, delta=0.0001)

        # ensure that both networks are isomorphic - same number of
        # nodes and same edge connection pattern.
        l_isomorphic = nx.is_isomorphic(mol_graph, test_molecule.topology)
        self.assertTrue(l_isomorphic)

        # if the network is isomorphic check that the nodes
        # are correctly labelled
        if l_isomorphic:
            for l_atom, l_test_atom in zip(mol_graph.nodes(),
                                           test_molecule.topology.nodes()):
                # check the atom types and ids are the same in both graphs
                self.assertEqual(l_test_atom.symbol,
                                 ['N', 'C', 'C', 'O'][l_test_atom.id])
                self.assertEqual(l_atom.symbol,
                                 ['N', 'C', 'C', 'O'][l_atom.id])

            # ensure the edges in the test molecule connect the same atoms
            # as the edges in the constructed molecule
            bond_labels = {tuple(sorted([bond[0].id, bond[1].id])) for bond in mol_graph.edges_iter()}
            test_bond_labels = {tuple(sorted([test_bond[0].id, test_bond[1].id])) for test_bond in test_molecule.topology.edges_iter()}

            # check that the sets are equivalent (a permutation of each other)
            self.assertTrue(test_bond_labels == bond_labels)

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

    def test_bond(self):
        bond = Bond(1, 2, 3, 4)
        self.assertEqual(1, bond.id)
        self.assertEqual(2, bond.atom1_id)
        self.assertEqual(3, bond.atom2_id)
        self.assertEqual(4, bond.order)


if __name__ == '__main__':
    unittest.main()
