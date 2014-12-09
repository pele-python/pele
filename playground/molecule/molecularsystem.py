from pele.systems.basesystem import BaseSystem
from playground.molecule.parser import Parser, PymolParser
from playground.molecule.molecule import Molecule, Protein

class MolecularSystem(BaseSystem):
    """ A container for general molecules which are defined as single component graphs.
    Maintains a dictionary of generic molecular objects and makes sure that
    each molecule has a unique Id within the collection which is the key in the
    dictionary.
    Also creates an instance of a parser which converts a file into a set of graphs.
    Each resulting graph is used to create a molecular instance.
    Molecules are defined as covalently linked atoms. So e.g separate protein chains
    count as distinct molecules.  
    """
    
    def __init__(self):
        # Initialise the member variables
        self.molecules = {}
        self.mol_id = 0 # every time a molecule is added increment this integer
        self.molecule_types = {}

    def select_parser(self, parser='Pymol', args='-qc'):
        # set up the parser for the system
        if parser == 'Pymol':
            self.parser = PymolParser(args)
        elif parser == 'Basic':
            self.parser = Parser()
        else:
            # Default parser is PymolParser
            self.parser = PymolParser(args)
        
    def load_file(self, l_filename, parser='pymol', moleculeType='general', args='-qc'):
        ''' Function looks at a file and adds every molecule in the file 
        to the molecular system as a separate molecule. Can specify the parser used
        to analyse the file.'''
        
        # create the specified parser object
        self.select_parser(parser, args)

        # parse the file to construct graphs of all the molecules in the file
        self.parser.parse(l_filename)
        
        # load all the coords and graphs obtained from this file into lists.
        l_graph_list = self.parser.get_graph_all()
        l_coords_list = self.parser.get_coords_all()
        
        # loop through the graphs returned from the parser
        # and add a molecule to the dictionary each time. 
        for graph, coords in zip(l_graph_list, l_coords_list):
            self.add_molecule(Molecule(self.mol_id, coords, graph))


    def add_molecule(self, molecule):
        ''' adds a molecule to the molecular system and gives it a system id.
        Also builds up a profile of how many of each type of molcules with distinct 
        hash_values there are.'''
        self.molecules[self.mol_id] = molecule
        molecule.change_id(self.mol_id)
        self.mol_id += 1
        
        # builds up a list of molcule ids that have a certain hash value
        try:
            self.molecule_types[molecule.hash_value] = [self.molecule_types[molecule.hash_value], molecule.id]
        except KeyError:
            self.molecule_types[molecule.hash_value] = [molecule.id]
        
        print self.molecule_types

    def __iter__(self):
        ''' Returns an iterator object over the molecules in the system.'''
        return self.molecules.itervalues()

    def __eq__(self, other):
        ''' Function determines the equality of two molecular systems
            based on the equality of the molecules within, and 
            the equality of the atoms with in those.
            Makes heavy use of the hash_value for a molecule which is expected to be unique 
            for most molecules.'''
        
        #assume both systems are equal
        equal = True 
        
        # get a list of the molecules in "self"
        l_mols_self = [ mol for mol in self.__iter__()]
        
        # attempt to get a list of the molecules "other"
        try:
            l_mols_other = [ mol for mol in other.__iter__()]
        except:
            # the iterator does not exist so is definitely not a molecular system.
            equal = False

        # don't do other checks if we already found non-equivalence            
        if equal:
            # check that the number of molecules in both lists is the same
            equal = (len(l_mols_self) == len(l_mols_other))

        if equal:
            # convert the list of molecules to sets - eliminates any molecules that have the same hashing function.
            # e.g. all copies of the same solvent molecule reduce to one entry in the set 
            mols_self = set(l_mols_self)
            mols_other = set(l_mols_self)
            
            # check that the types of molecule in the two systems are the same. i.e. that neither set has a 
            # molecule which is not in the other set. 
            # To help with this the hashing function in the molecule class has been over ridden.
            # The hash function for each molecule object is now the product of the absolute values of the non-zero eigenvalues 
            # of the adjacency matrix where the edges are weighted by the product of atomic weights at either end of the bond.
            # Identical molecules will generate the same hash. There may be one or two random esoteric cases that
            # generate the same hash, but this will be fluke, and the graphs are unlikely to be isomorphic.
            # Comparing the sets checks that the set of hash_functions is the same, and if they are the same
            # then the element equivalence operator (molecule.__eq__ function) is called which checks that the molecules 
            # have isomorphic graphs. Nodes at equivalent points in the isomorphics graph are compared using the node attribute "atom"
            # This is a reference to the Atom objects associated with each node which triggers the Atom.__eq__ function 
            # which checks that equivalent nodes have the same label. (a lot going on behind one simple statement!!)
            equal = mols_self == mols_other
        
        if equal:
            # we already know that the distinct number of types of molecules and the numerical values of the hash functions
            # of those classes of molecule are the same for both systems. Now we check the number of each type of molecule.
            # Which we count as we add/delete molecules from the system. 
            # Convert one of the sets of molecules to a list and extract the hash value of each unique type of molecule in the system.
            # Check that the number of each type of molecule is the same in both systems. 
            # Probably more relevant for solvent molecules and that sort of thing.
            if False in [len(self.molecule_types[h]) == len(other.molecule_types[h]) for h in [mol.hash_value for mol in list(mols_self)]]:
                equal = False
        
        return  equal

class ProteinSystem(MolecularSystem):
    '''  A class for handling systems of pure proteins. For ligands and other
         kinds of mixed systems create your own class for dealing with that.
         essentially this is intended to handle multi chain proteins. '''
    def __init__(self):
        super(ProteinSystem, self).__init__()

    def load_file(self, l_filename, parser='pymol', moleculeType='general', args='-qc'):
        ''' Function overloads the function of the same name in the base class.
        It looks at a file and adds every molecule in the file 
        to the molecular system as a separate molecule, assuming the molecules 
        are all proteins. Can specify the parser used
        to analyse the file. The Protein class
        checks for the basic nature of the protein and raises an exception 
        if a molecule is loaded that is not a peptide polymer.'''
        
        # create the specified parser object
        self.select_parser(parser, args)

        # parse the file to construct graphs of all the molecules in the file
        self.parser.parse(l_filename)
        
        # load all the coords and graphs obtained from this file into lists.
        l_graph_list = self.parser.get_graph_all()
        l_coords_list = self.parser.get_coords_all()
        
        # loop through the graphs returned from the parser.
        for graph, coords in zip(l_graph_list, l_coords_list):
            self.molecules[self.mol_id] = Protein(self.mol_id, coords, graph)
            self.mol_id += 1
    
