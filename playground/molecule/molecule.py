import pele.utils.elements as elem

import networkx as nx
import pymol_wrapper
import os as os


class Molecule(object):
    """
    A representation of a molecule consisting of the following items:

    1) A graph whose nodes are atom objects and whose edges have bonds
       attributes. This defines the topology of the molecule.
    2) A list of coords of the atomic positions defining the current
       state of the molecule.
    3) A collection of functions for performing operations on and
       calculations about the molecule

    self.topology     Network X Graph. Nodes are atoms and edges are bonds.
    self.coords       Spatial coords of the atoms. The unique index number, i,
                      of each node in the network corresponds to the ith
                      entry in the coords array.
    """

    def __init__(self, l_input_filename, args):
        '''
        Initialises an instance of a molecular class from a file of data.

        Takes advantage of the the pymol file interface to parse a range of
        files such as PDBs and XYZ files

        Constructs a graph of the molecule.
        and maintains a list of atomic coords.
        '''
        # Set up Pymol
        pymol = pymol_wrapper.init(args)

        # load the data into pymol
        self.load_data_pymol(pymol, l_input_filename)

        # Extract the pymol object defining the molecule
        l_pymol_data = self.extract_data_pymol(pymol)

        # Initialise member variables of Molecule object
        self.coords = []
        self.topology = nx.Graph()

        # Gget the atomic information from the pymol data.
        # This populates the self.coords array and
        # constructs the nodes of self.topology but with no edge information.
        l_pymol_index_map = self.extract_atoms(l_pymol_data)

        # Construct edges of the graph from the pymol data and the index map.
        # Populates self.topology.
        self.extract_topology(l_pymol_data, l_pymol_index_map)

    def load_data_pymol(self, pymol, l_input_filename):
        '''Loads a file into pymol object to take advantage of it's file handling
        and molecular construction capability.'''

        try:
            with open(l_input_filename, "r") as f_in:
                f_in.close()
        except IOError:
            raise Exception("Unable to open file: " + l_input_filename)

        try:
            # Determine file type from file header
            l_fileType = l_input_filename[-3:]
        except IndexError:
            raise Exception("Cannot find file header: " + l_input_filename)

        # Process pdb file
        if l_fileType == 'pdb':
            # Attempt to load the data
            pymol.cmd.load(l_input_filename)
        elif l_fileType == 'xyz':
            # Pymol loads all the frames in an xyz file,
            # This function loads only the first frame.
            try:
                with open(l_input_filename, "r") as f_in:
                    l_raw_data = f_in.readlines()
            except IOError:
                raise Exception("Cannot open input file:" + l_input_filename)

            # Extract the number of atoms
            l_num_atoms = int(l_raw_data[0])

            try:
                with open('temporaryFile.xyz', "w") as f_out:
                    f_out.writelines(l_raw_data[0:l_num_atoms + 2])
            except IOError:
                raise Exception("Unable to write xyz frame for output")

            pymol.cmd.load('temporaryFile.xyz')
            os.remove('temporaryFile.xyz')
        else:
            raise Exception("Unrecognised file type: " + l_input_filename)

    def extract_data_pymol(self, pymol):
        '''Function extracts all the data about a single model from pymol'''
        # Get the name of the last object loaded into pymol instance
        # and return all information about that pymol object
        try:
            pymol_object_name = pymol.cmd.get_object_list()[-1]
        except IndexError:
            raise Exception("No objects in list obtained from pymol")
        return pymol.cmd.get_model(pymol_object_name)

    def extract_atoms(self, l_pymol_object_data):
        '''Extract the coordinates of each atom from pymol
           Constructs an atom object for each atom and adds it as a node
           to the graph. Returns a dictionary which maps the pymol index to the
           pele index.'''

        # Initialise index to keep track of position in coords array
        l_atom_id = 0

        # Initialise a dictionary mapping the pymol id to the pele atoms
        l_pymol_index_map = {}

        # Loop through the atoms in the pymol object and convert into
        # The pele molecular representation.
        for l_atom_pymol in l_pymol_object_data.atom:
            try:
                # Populate the coords array coords
                self.coords.append(l_atom_pymol.coord[0])
                self.coords.append(l_atom_pymol.coord[1])
                self.coords.append(l_atom_pymol.coord[2])
            except LookupError:
                raise Exception("Error extracting data from Pymol Object")

            # Create a pele atom object
            l_atom = Atom(l_atom_id, l_atom_pymol.name)

            # Add the current atom as a node
            self.topology.add_node(l_atom)

            try:
                # Make a note of the reverse mapping between the the pymol
                # index and the pele atom.
                l_pymol_index_map[l_pymol_object_data.index_atom(l_atom_pymol)] = l_atom
            except LookupError:
                raise Exception("Invalid key for atom dictionary")

            # Increment atomic ID.
            l_atom_id += 1

        return l_pymol_index_map

    def extract_topology(self, l_pymol_object_data, l_pymol_index_map):
        '''This function loops through the bond information from pymol,
           identifies the pele atoms referred to in the pymol bond and adds
           the edge to the graph. A pele bond object is associate with the
           edge. '''
        # Initialise index to keep track of bond id
        l_bond_id = 0

        # Loop through the pymol bond data
        for l_bond in l_pymol_object_data.bond:
            try:
                # Find the pele atom referred to in the pymol bond info
                l_atom1 = l_pymol_index_map.get(l_bond.index[0], None)
                l_atom2 = l_pymol_index_map.get(l_bond.index[1], None)
            except LookupError:
                raise Exception("Unable to find Pele atom with Pymol index.")

            # Make sure that both atoms were found
            if None not in [l_atom1, l_atom2]:
                # Create a pele bond object
                bond = Bond(l_bond_id,
                            l_atom1.id,
                            l_atom2.id,
                            l_bond.order)

                # Add the edge to the network.
                # Associate the bond object to the edge.
                self.topology.add_edge(l_atom1, l_atom2, object=bond)

                # Increase count
                l_bond_id += 1


class Atom(object):
    '''
    A minimal representation of an atom:
        id is unique index in the molecule.
        symbol identifies the type of atom.

        Conventions for atomic symbols that are understood are the PDB
        atom naming system and the normal element naming system.
    '''
    def __init__(self, id, symbol):
        # populate the member variables of the Atom class.
        self.id = id

        # A variety of atomic naming conventions within molecules.
        # exist. E.g. H could be 1H, 2H, OH etc in a pdb file.
        # Keep a note of the given defining symbol.
        self.alt_symbol = symbol

        # lookup the standard symbol in the alternative names dictionary.
        try:
            self.symbol = elem.alternate_names[self.alt_symbol]
        except KeyError:
            # if there is no entry in the alt_names list for the identifier,
            # assume the identifier is the alternative_symbol.
            self.symbol = self.alt_symbol

        # look up the remaining element details from the elements dictionary
        try:
            self.name = elem.elements[self.symbol]['name']
        except KeyError:
            # Raise an exception pointing out that the given symbol is unknown.
            raise Exception("Unknown element: " + symbol)

        # given symbol has been identified. Acquire remaining details
        self.mass = elem.elements[self.symbol]['mass']
        self.radius = elem.elements[self.symbol]['radius']
        self.color = elem.elements[self.symbol]['color']


class Bond(object):
    '''
    defines a bond between two atoms within a molecule. Atom1 and atom2
    are the unique ids of two atoms within a molecule and bond_id is the
    unique id of a given bond.
    '''
    def __init__(self, id, atom1_id, atom2_id, order):
        self.id = id
        self.atom1_id = atom1_id
        self.atom2_id = atom2_id
        self.order = order

if __name__ == '__main__':

    HETS = Molecule('2RNM HETS.pdb', '-qc')
    PHG12 = Molecule('lowest.xyz', '-qc')

    print "Done"
