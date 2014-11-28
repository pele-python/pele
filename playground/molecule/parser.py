from playground.molecule.molecule import Atom


import networkx as nx
import numpy as np
import abc


class Parser(object):
    '''The parser class and its subclasses take 
    a file name containing information
    about a system of one or more molecules and 
    produces a graph whose nodes are atom objects and whose
    edges represent the covalent bonds between the atoms.

    Also the initial coords of the atoms are stored in an array
    and indexed by the atom.id.
       
    Each frame or model in the input file is analysed for connectivity 
    separately, but is added as a disconnected components to the same graph.
    
    The entire graph of all frames can then be analysed to find 
    all the parts of the network which are not bonded to each other.

    This automatically takes care of multiple frames and models, 
    multiple chains within each model as well as solvents 
    molecules, hetatoms, ligands etc.  If groups of atoms in the 
    molecular network are not covalently linked then they become a distinct 
    regions of the graph and ultimately distinct molecules.  

    After the main graph has been split the coords array for each sub_graph
    is reconstructed. The atoms within each sub graph are relabelled so
    atom.id indexes their coords in their own coords list. 
    '''
    __metaclass__ = abc.ABCMeta

    def __init__(self, args):
        self.args = args
        self.topology = nx.Graph()
        self.coords = []
        
    def parse(self,filename):
        ''' FUnction manages the top level process of parsing a file'''
        self.input_filename = filename
        self.load_data()
        self.form_graph()
        self.find_distinct_graphs()

    def load_data(self):
        ''' use the filename variable and args to load data.'''
        try:
            with open(self.input_filename, "r") as f_in:
                self.raw_data = f_in.readlines()
        except IOError:
            raise Exception("Unable to open file: " + self.input_filename)

    def form_graph(self):
        ''' analyse the raw data to construct a graph.
            Each frame or model must be analysed separately.
            Each node in the graph must consists of an Atom object.
            A unique id must be assigned to each atom which corresponds to the coords data index.
            The coords is the standard [x1, y1, z1, x2, y2, z2,...] arrangement. '''
        return

    def get_graph_all(self):
        ''' Access function returns a list of molecular graphs whose 
        nodes are Atom objects and edges are bonds.'''
        return self.component_list
     
    def get_coords_all(self):
        ''' Access function returns the coords of the nodes as an array.
        Each Atom class has a unique id which doubles as the index of
        the coords in the coords array. 
         '''
        return self.coords_list

    def find_distinct_graphs(self):
        ''' Analyse the components of a graph to find all the disconnected components.
        Each list of nodes is converted into a sub graph and the coords information is broken up 
        into corresponding sublists.'''
 
        # analyse the components of the topology to find all the disconnected components
        self.component_list = nx.connected_component_subgraphs(self.topology)
 
        # initialise a list of coords lists to have the same number of elements as there are components
        self.coords_list = [None] * len(self.component_list)
 
        # loop through the components
        for l_component_index, l_component in enumerate(self.component_list):
 
            # construct a sublist with three times the number of entries as there are atoms in the component
            self.coords_list[l_component_index] = [None] * 3 * len(l_component)
 
            # loop through the atoms in the component and extract its coords from the original coords array,
            # which is indexed by the atoms original id.
            # place it in the new sublist at the new index and renumber the atom with the new index.
            l_new_atom_index = 0
            for atom in l_component.nodes():
                # transfer the current atoms coords to the right place in the new array
                self.coords_list[l_component_index][3 * l_new_atom_index] = self.coords[3 * atom.id]
                self.coords_list[l_component_index][3 * l_new_atom_index + 1] = self.coords[3 * atom.id + 1]
                self.coords_list[l_component_index][3 * l_new_atom_index + 2] = self.coords[3 * atom.id + 2]
 
                # renumber atom 
                atom.change_id(l_new_atom_index)
 
                # increment the new atom index 
                l_new_atom_index += 1
 
class PymolParser(Parser):
    ''' This class takes an xyz or pdb file, and uses the pymol interface to
        parse the file and construct bond information.
        
        A graph is then constructed from the pymol information with
        a list of coords for each atom whose order correspondings to the
        id number in the atom class of each node. 
        '''

    def __init__(self, args='-qc'):
        ''' Initialise the pymol parser.'''
        # initialise the base class
        super(PymolParser, self).__init__(args)

        # Load the pymol instance '''
        try:
            import __main__

            __main__.pymol_argv = ['pymol', args]
            import pymol

            self.pymol = pymol
            self.pymol.finish_launching()
        except:
            raise Exception("Unable to launch pymol")
        
    def parse(self, filename):
        ''' PymolParser overload the Base Object Parse routine.'''
        self.filename = filename
        # load the data into pymol
        self.load_data()
        # Extract the pymol object defining the molecule
        self.extract_data_pymol()
        # form the main graph of all the data in the file
        self.form_graph()
        # split the main graph into separate molecules
        self.find_distinct_graphs()

    def load_data(self):
        '''Loads a file into pymol to take advantage of it's file handling
        and molecular construction capability.'''
        try:
            with open(self.filename, "r") as f_in:
                f_in.close()
        except IOError:
            raise Exception("Unable to open file: " + self.filename)

        try:
            # Determine file type from file header
            self.fileType = self.filename[-3:]
        except IndexError:
            raise Exception("Cannot find file header: " + self.filename)

        # Process file
        if self.fileType in ['pdb', 'xyz']:
            try:
                # Attempt to load the data
                self.pymol.cmd.load(self.filename)
            except:
                raise Exception("Pymol unable to load filename: " + self.filename)
        else:
            raise Exception("Unrecognised filetype: " + self.fileType)

    def extract_data_pymol(self):
        '''Function extracts all the data about a single model from pymol'''
        # Get the name of the last object loaded into pymol instance
        # and return all information about that pymol object
        try:
            self.pymol_object_name = self.pymol.cmd.get_object_list()[-1]
        except IndexError:
            raise Exception("No objects in list obtained from pymol")

        # extract the data and load it into a member variable
        self.pymol_data = self.pymol.cmd.get_model(self.pymol_object_name)

    def form_graph(self):
        ''' Function extracts the data from the pymol object 
            and uses it to construct a graph in two stages.
            
            First, extract the coordinates of each atom from pymol
            create a unique id for the atom 
            and use it to construct an atom object for each atom. 
            Add the object as a node to the main graph. 
            Create a dictionary index by the pymol id that returns 
            a reference to the atom object.
            
            Second, extract the bond information from the pymol
            and use the dictionary to determine 
            which pele atoms the bond refers to and use 
            this to construct an edge between the two atoms.'''

        # initialise the pymol index map
        l_pymol_index_map = {}

        # Initialise index to keep track of position in coords array
        l_atom_id = 0

        # Loop through the atoms in the pymol object and convert into
        # The pele molecular representation.
        for l_atom_pymol in self.pymol_data.atom:
            try:
                # Populate the coords array coords
                self.coords.append(l_atom_pymol.coord[0])
                self.coords.append(l_atom_pymol.coord[1])
                self.coords.append(l_atom_pymol.coord[2])
            except LookupError:
                raise Exception("Error extracting data from Pymol Object")

            # Create a pele atom object
            l_atom = Atom(l_atom_id, l_atom_pymol.name)

            # Add the current atom as a node, with the node attribute symbol
            self.topology.add_node(l_atom, atom=l_atom)

            try:
                # Make a note of the reverse mapping between the the pymol
                # index and the pele atom.
                l_pymol_index_map[self.pymol_data.index_atom(l_atom_pymol)] = l_atom
            except LookupError:
                raise Exception("Invalid key for atom dictionary")

            # Increment atomic ID.
            l_atom_id += 1

        # Loop through the pymol bond data
        for l_bond in self.pymol_data.bond:
            try:
                # Find the pele atom referred to in the pymol bond info
                l_atom1 = l_pymol_index_map.get(l_bond.index[0], None)
                l_atom2 = l_pymol_index_map.get(l_bond.index[1], None)
            except LookupError:
                raise Exception("Unable to find Pele atom with Pymol index.")

            # Make sure that both atoms were found and add edge to the network
            if None not in [l_atom1, l_atom2]:
                # Add the edge to the network 
                self.topology.add_edge(l_atom1, l_atom2, hweight=np.sqrt(l_atom1.mass * l_atom2.mass))

