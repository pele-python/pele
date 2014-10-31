from utilities import read_text_file, write_text_file

import networkx as nx
import numpy as np
import pymol_wrapper
import sys as sys
import os as os

class Molecule():
    """
    A representation of a molecule consisting of the following items:
    
    1) A graph whose nodes are atom objects and whose edges have bonds attributes.
        This defines the topology of the molecule.
    2) A list of coords of the atomic positions defining the current state of the molecule.
    3) A collection of functions for performing operations on and calculations about the molecule
    
    
    self.topology              Network X Graph. The nodes are atoms and the edges are bonds.
    self.coords                Spatial coords of the atoms. The unique index number, i, of each node in 
                               the network corresponds to the ith entry in the coords array.  
    
    """

    def __init__(self, l_input_filename):
        '''
        Initialises an instance of a molecular class from a file of data.
        
        Takes advantage of the the pymol file interface to parse a range of
        files such as PDBs and XYZ files
        
        Constructs a graph of the molecule.
        and maintains a list of atomic coords            
        '''
        #Set up Pymol
        pymol=pymol_wrapper.init('-qc')
            
        #load the data into pymol    
        self.load_data_pymol(pymol,l_input_filename)

        #Extract the pymol object defining the molecule
        l_pymol_data=self.extract_data_pymol(pymol)
        
        #initialise member variables of Molecule object
        self.coords=[]
        self.topology=nx.Graph()
        
        #get the atomic information from the pymol data.
        #This populates the self.coords array and 
        #constructs the nodes of self.topology but with no edge information.
        l_pymol_index_map=self.extract_atoms(l_pymol_data)
        
        #construct the edges of the graph from the pymol data and the index map.
        #Populates self.topology 
        self.extract_topology(l_pymol_data,l_pymol_index_map)
        

    def load_data_pymol(self,pymol,l_input_filename):
        '''Loads a file into pymol object to take advantage of it's file handling
        and molecular construction capability.'''
        try:
            #determine file type from file header
            l_fileType=l_input_filename[-3:]
            
            #process pdb file
            if l_fileType=='pdb':
                #attempt to load the data
                pymol.cmd.load(l_input_filename)

            if l_fileType=='xyz':
                #pymol loads all the frames in an xyz file,
                #This function loads only the first frame.
                l_raw_data=read_text_file(l_input_filename)
                l_num_atoms=int(l_raw_data[0])
                write_text_file('temporaryFile.xyz',l_raw_data[0:l_num_atoms+2-1])
                pymol.cmd.load('temporaryFile.xyz')
                os.remove('temporaryFile.xyz')

        except: 
            #unable to load file
            print "Unable to load filename: ",l_input_filename
            sys.exit()
    
    def extract_data_pymol(self,pymol):
        '''Function extracts all the data about a single model from pymol''' 
        try:
            #Get the name of the last object that was loaded into the local pymol instance
            #and return all information about that pymol object
            return(pymol.cmd.get_model(pymol.cmd.get_object_list()[-1]))  
        except:
            print "Unable to extract molecular information from Pymol"
            sys.exit()
    
    def extract_atoms(self,l_pymol_object_data):
        '''Extract the coordinates of each atom from pymol 
           Constructs an atom object for each atom and adds it as a node to the graph.
           returns a dictionary which maps the pymol index to the pele index.'''
        
        try:
            #index keeping track of position in coords array
            l_atom_id=0 
            l_pymol_index_map={}
            #loop through the atoms in the pymol object and convert into 
            #the pele molecular representation.
            for l_atom_pymol in l_pymol_object_data.atom:
                    #populate the coords array coords
                    self.coords.append(l_atom_pymol.coord[0])
                    self.coords.append(l_atom_pymol.coord[1])
                    self.coords.append(l_atom_pymol.coord[2])
                    
                    #create a pele atom object
                    l_atom=Atom(l_atom_pymol.name, l_atom_id)
                    
                    #Add the current atom as a node
                    self.topology.add_node(l_atom)
                    
                    #make a note of the reverse mapping between the the pymol index and the pele atom.
                    l_pymol_index_map[l_atom_pymol.id]=l_atom
                    
                    #increment atomic ID.
                    l_atom_id+=1
                    
            return l_pymol_index_map
        except:
            print "Unable to process atomic information from pymol"
            sys.exit()

    def extract_topology(self,l_pymol_object_data, l_pymol_index_map):
        '''This function loops through the bond information from pymol, identifies 
           the pele atoms referred to in the bond and adds the edge to the graph.
           Creates a pele bond object and associates it with the edge.
        '''
        
        #loop through the bond data
        try:
            l_bond_id=0
            for l_bond in l_pymol_object_data.bond:
                    #find the pele atom referred to in the pymol bond info 
                    l_atom1=l_pymol_index_map.get(l_bond.index[0],None)
                    l_atom2=l_pymol_index_map.get(l_bond.index[1],None)

                    if not None in [l_atom1,l_atom2]: #for some reason pymol reports a bond between non-existent atoms...??
                        #create a pele bond object
                        bond=Bond(l_bond_id,l_atom1.atom_id,l_atom2.atom_id,l_bond.order)
                        
                        #add the edge to the network and associate the bond object to the edge
                        self.topology.add_edge(l_atom1,l_atom2,object=bond)

                        #increase count
                        l_bond_id+=1
                    
        except:
            print "Unable to process bond information from pymol"
            sys.exit()

        print "topology done"

class Atom(object):
    '''
    A minimal representation of an atom:
        Mass, atomic number, name, symbol
         
    '''
    def __init__(self, atom_name, atom_id):
        self.atom_name=atom_name
        self.atom_id=atom_id
        #self.atom_symbol=atom_symbol
        #self.atom_mass=elem[self.atom_symbol]['mass']
        #self.atom_number=elem[self.atom_symbol]['atomicNumber']

class Bond(object):
    ''' 
    defines a bond between two atoms within a molecule. Atom1 and atom2 are the unique ids 
    of two atoms within a molecule and bond_id is the unique id of a given bond
    '''
    def __init__(self,bond_id,atom1_id,atom2_id,bond_order):
        self.bond_id=bond_id
        self.atom1_id=atom1_id
        self.atom2_id=atom2_id
        self.order=bond_order

        
if __name__ == '__main__':
    
    HETS=Molecule('2RNM HETS.pdb')
    PHG12=Molecule('lowest.xyz')

    print "Done"
    
            