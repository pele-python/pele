import pele.utils.elements as elem

from operator import mul
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
import os as os




class Atom(object):
    """
    A minimal representation of an atom:
        id is unique index in the molecule.
        symbol identifies the type of atom.

        Conventions for atomic symbols that are understood are the PDB
        atom naming system and the normal element naming system.
    """

    def __init__(self, atom_id, symbol):
        # populate the member variables of the Atom class.
        self.id = atom_id

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

    def change_id(self, new_id):
        self.id = new_id

    def __eq__(self, other):
        ''' Function to determine if two atoms are the same type'''
        try:
            equal = self.symbol == other.symbol 
        except:
            # if other does not have the attribute symbol then it is not an Atom class
            # so cannot be equal!
            equal = False
        
        return equal

    def __ne__(self, other):
        ''' Function to determine if two atoms are not same type'''
        return self.symbol != other.symbol

class Molecule(object):
    """
    A representation of a molecule consisting of the following items:

    1) A graph whose nodes are atom objects and whose edges have bonds
       attributes. This defines the topology of the molecule.
    2) A list of coords of the atomic positions defining the current
       state of the molecule.
    3) A collection of functions for performing operations to 
       molecules, such as rotate, translate etc.


    self.topology     Network X Graph. Nodes are atoms and edges are bonds.
    self.coords       Spatial coords of the atoms. The unique index number, i,
                      of each node in the network corresponds to the 
                      position of its coords in the array.
                      So atom.id x coord is as coords[3*atom.id], 
                                          y at coords[3*atom.id + 1]
                                          z at coords[3*atom.id + 2]
    """

    def __init__(self, mol_id, coords, topology):
        """ Initialises an instance of a molecular class """
        self.id = mol_id
        self.topology = topology
        self.coords = coords
        
        # compute the hash value at initialisation
        self.hash_value = None
        self.__hash__()

    def change_id(self, new_id):
        ''' function re identifies the molecule '''
        self.id = new_id

    def __hash__(self):
        ''' Define a hashing function for the molecule.
        Computes the product of the non-zero eigenvalues of a weighted adjacency matrix.
        The weights are the square root of the products of the atomic weights
        of the connected vertices (computed at network creation).'''
        # if hash_value has already been computed then don't compute it again.
        if self.hash_value == None:
            #compute the eigenvalues of the adjacency matrix weighted by hweight.
            evals = np.linalg.eigvals(nx.adjacency_matrix(self.topology, weight='hweight'))
            
            # compute and store the product of absolute values of eigenvalues that are higher than delta_e
            self.hash_value = reduce(mul, [abs(e) for e in evals if abs(e)>0.001] , 1)
            
            # truncate the value to 3 decimal place
            self.hash_value = float(int(self.hash_value * 1000))/1000.0
            
        return self.hash_value
    
    def __eq__(self, other):

        ''' Function has two molecules as input and decides whether or not 
            they are the same molecule based on the graph. 
            Both the node labels and graph topology must match '''

        # set up the function that is called to determine the node attribute
        # used to determine equality of the nodes        
        nm = iso.categorical_node_match('atom', Atom(0, 'Xx'))
        
        # check that both networks are isomorphic - same number of
        # nodes and same edge connection pattern and that the node atom types 
        # are equivalent 
        return nx.is_isomorphic(self.topology, other.topology, node_match=nm)

    def __ne__(self, other):

        ''' Function has two molecules as input and decides whether or not 
            they are the same molecule based on the graph. 
            Both the node labels and graph topology must match to be equal '''

        # set up the function that is called to determine the node attribute
        # used to determine equality of the nodes        
        nm = iso.categorical_node_match('symbol','H')
        
        # check that both networks are isomorphic - same number of
        # nodes and same edge connection pattern and that the node atom types 
        # are equivalent 
        return nx.is_isomorphic(self.topology, other.topology, node_match=nm)
    

class Protein(Molecule):
    """ A sub class of the Molecule base class.
        which contains objects that relate to solely to proteins.
        The idea is that the Protein Class contains only a single 
        covalently linked chain.  This is because it is a sub
        class of Molecule which is designed to represent only a single chain.
        For multiple chains use the ProteinSystem Class.

        self.name - a unique name for the protein.
        self.c_alphas - a list of the c_alphas atoms for each chain
        self.residues  - A list of residue objects for each chain
        self.peptide_bonds - a list of sub graphs representing each peptide bond in each chains
        self.sidechains - a list of sub graphs of all the side chains in each chain
        self.backbone  - the correct entire backbone for each chain excluding termini
        self.mainchain  - the entire main chain of each chain

        self.dihedrals:  The backbone dihedral angles are computed as a list of dynamically changing information.

        """

    def __init__(self, protein_id, coords, topology):
        # initialise the base class

        super(Molecule, self).__init__(protein_id, coords, topology)

        # Identify the longest path in each chain. This is mostly the backbone, 
        # but may include terminal groups or side chain of the penultimate residue
        # depending on which gives the longest path.
        self.find_longest_paths()

        # Identify the number of cliques to which each atom belongs in the molecule.
        # Store grouped as chains
        self.find_cliques()

        # create a list of four clique carbon atoms in each of the longest chains. 
        # this may includes carbons in any terminating group or
        # side chains of the terminating residues.
        self.find_four_clique_carbons()

        # analyse the neighbourhood of each the four_clique_carbons to determine 
        # sub graphs for each terminus and residue. 
        # Assign a unique id to each residue.
        # within each residue determine the sub graphs for the 
        # back bone atoms, main _chain and side_chain.
        self.generate_residues_from_four_clique_carbons()

    def generate_residues_from_four_clique_carbons(self):
        """ Takes each carbon alpha in the carbon alpha list and analyses it to see if it belongs
            to a residue, a terminating group, or is part of the side chain of the terminating residue.
            A residue object is then created with a unique id for each residue in the protein.
            The residue object has sub graphs of the entire residue, the side_chain,
            all the backbone atoms, and  main chain atoms.
            Terminating groups are treated as residues in their own right with their own IDs
            but contain only backbone atoms or main_chain groups.

            At the end of the process if there are any atoms in the molecule which have not been
            allocated to a standard residue object an exception is raised.
        """

        # check there are the same numbers of chains as there are lists of carbon_alphas. 
        if len(self.chains) != len(self.four_clique_carbons):
            raise Exception("Number of lists of four clique carbons doesn't match number of chains")

        # set the first res_id
        res_id = 1

        # loop through each chain and list of c_alphas
        for chain_index, four_clique_carbons in enumerate(self.four_clique_carbons):
            # loop through the four_clique_carbons in this chain
            for four_clique_carbon in four_clique_carbons:
                # analyse the neighbour hood of the four clique carbon.
                outgraph = self.analyse_four_clique_carbon(four_clique_carbon, chain_index)

    def draw_local_neighbourhood(self, atom, depth, graph_type='spring'):
        """ Function draws the neighbours of a graph around a node,
            and then the neighbours of those neighbours.  The depth of the drawing defines
            how many generations of neighbours are plotted. Each generation is
            given a distinct colour."""
        # start the output list with the seed atom.
        output_list = [atom]

        # generate a dictionary to store the node colour map for the network.
        node_colour_dict = {}

        # store the root atom colour as green. Use the atom id as a dictionary key.
        node_colour_dict[atom.id] = [0.0, 1.0, 0.0]

        # loop through the depths (i.e how many neighbours away we visit from the root node).
        for cur_depth in range(0, depth + 1):

            # create a temporary copy of the output list so far
            cur_list = output_list[:]  # make a copy of output list

            # loop through the entire temporary output list
            for node in cur_list:

                # for each node in the output list get all of it's neighbours from the main topology graph. 
                neighbours = nx.all_neighbors(self.topology, node)

                # loop through each neighbour in the current output list.
                for n in neighbours:

                    # if the current atom is not anywhere in the output list already then it is accepted.  

                    if not n in output_list:
                        output_list.append(n)

                        # a colour is assigned depending on which depth level we have achieved.
                        if cur_depth == 0:
                            node_colour_dict[n.id] = [1.0, 0.0,
                                                      1.0]  # magenta for the first level node (first time through!)
                        if cur_depth == 1:
                            node_colour_dict[n.id] = [1.0, 1.0, 0.0]  # yellow for the second level node
                        if cur_depth == 2:
                            node_colour_dict[n.id] = [0.0, 1.0, 1.0]  # cyan for the third level nodes
                        if cur_depth > 2:  # increasing greyness from white to black over 10 levels.
                            node_colour_dict[n.id] = [1.0 * float(13 - cur_depth) / float(10.0),
                                                      1.0 * float(13 - cur_depth) / float(10.0),
                                                      1.0 * float(13 - cur_depth) / float(10.0)]

        # define the sub graph induced on self.topology by output_list nodes.
        G = nx.subgraph(self.topology, output_list)

        # figure out the colour of the nodes in the new graph - 
        # subgraph creates a copy of the nodes so using the unique id generated by the molecule class 
        # as the key for identify the nodes rather than the node itself. 
        node_colour_list = [node_colour_dict[n.id] for n in G.nodes_iter()]

        # plot the graph using the node_colours.
        self.draw_graph(G, graph_type=graph_type, node_size=150, node_colour=node_colour_list)

    def analyse_four_clique_carbon(self, four_clique_carbon, chain_index):
        """ The function takes a single four clique carbon and explores its neighbourhood.
            it returns a list of four graph objects:
                A residue_graph,
                sidechain_graph,
                backbone_graph,
                main_chain_graph

            if these objects cannot be determined then the corresponding entry in the list is None

        """

        # construct the output array

        # set up lists to hold the atoms in the various bits as we identify them
        back_bone_list = []
        main_chain_list = []
        side_chain_list = []
        residue_list = []

        # find the neighbourhood of C. Returns an iterator for the dictionary which we make into a list
        neighbours = [atom for atom in nx.all_neighbors(self.topology, four_clique_carbon)]

        # set up variables to store the backbone atoms as we encounter them
        H = None
        C = None
        O = None
        NH = None
        N = None
        CB = None
        HB = None

        self.draw_local_neighbourhood(four_clique_carbon, 3)
        print four_clique_carbon.id

        # loop through the four neighbours of the four clique carbon and perform four tests based on the 
        # neighbours
        for curr_neighbour in neighbours:

            # make list of the neighbour types of the current neighbour
            n_neighbours = [atom for atom in nx.all_neighbors(self.topology, curr_neighbour)]

            # Check the 4 clique carbon against the necessary conditions to be a carbon alpha.
            # Need to identify N, NH, C, O and CH, CB or HB
            # A carbon alpha will have:
            # 1) N as neighbour where the N is also a member of the longest chain and has at least one H as a neighbour
            # (NH)
            if curr_neighbour.symbol == 'N' and curr_neighbour in self.longest_paths[chain_index]:
                # Make a list of Hydrogen neighbours of the N in the main chain.
                NH_list = [nn for nn in n_neighbours if (nn.symbol == 'H')]

                if len(NH_list) == 2:
                    # if N is in the main chain and has two hydrogens then this residue could be the terminus.
                    NTerminus_flag = 1
                    N = curr_neighbour
                    NH = NH[0]
                    NH1 = NH[1]

                elif len(NH) == 1:
                    # Only one H and in the main chain suggests this is the NH of the peptide bond in a regular residues
                    NH = NH[0]
                    N = curr_neighbour
                elif len(NH) == 0:
                    # an N in the main chain with no H's
                    # check for case where all three of NN are carbons: this occurs in Proline!
                    NH = None

            # 2) C as neighbour where the C is also a member of the longest chain with a clique of 3 and one of it's
            # neighbours is an Oxygen.
            if curr_neighbour.symbol == 'C' and curr_neighbour in self.longest_paths[chain_index] and \
               self.cliques[chain_index][curr_neighbour] == 3:
                O = [nn for nn in n_neighbours if (nn.symbol == 'O')]
                if len(O) > 0:
                    O = O[0]
                    C = curr_neighbour
                else:
                    O = None

            # 3) X as a neighbour
            # where
            # X is either:
            # H in Glycine
            if curr_neighbour.symbol == 'H' and H is not None:  # i.e. H is already set and we have encountered a
            # second H.
                HB = curr_neighbour
            # or C in every other (standard) amino acid side chain (CB) which is always clique of 4 too.
            if curr_neighbour.symbol == 'C' and not curr_neighbour in self.longest_paths[chain_index] and \
               self.cliques[chain_index][curr_neighbour] == 4:
                CB = curr_neighbour

            # 4) H as a neighbour
            if curr_neighbour.symbol == 'H' and H is None:
                H = curr_neighbour

        # check to see if all the atoms in the back bone could be determined.
        if H and C and N and NH and O and (CB or HB):
            # now fairly sure that C is a carbon alpha in a protein. There may be other molecules with this arrangement
            # of stuff in other contexts, but we are allowed to assume the molecule is always a protein, in which case
            # this is very likely to be ok because I could find no other four clique carbons in a protein (with standard
            # amino acids) that matches these conditions exactly. Actually the following is a good test: will generate a
            # protein with at least one of every amino acid in it to make sure this is true!

            # The back bone is the N, C, NH and H and O identified above.
            # In a terminal amino acid sometimes this group will have an NH2 and sometimes an extra O and H.
            back_bone_list = [N, NH, four_clique_carbon, H, C, O]

            # The main chain is the N, C and target carbon identified above. 
            main_chain_list = [N, four_clique_carbon, C]

            # Determine list of nodes in side chains.
            # two special cases:  Prolines and Glycines.
            # Prolines are rings.
            # Glycines contain only HB
            # HB is not defined for any amino acid except glycine.

            if HB:
                side_chain_list = [HB]
            else:
                # Find the number of shortest paths between four_clique_carbon and CB.
                paths = [path for path in nx.all_shortest_paths(self.topology, four_clique_carbon, CB)]
                if len(paths) == 2:
                    # probably a proline.
                    print "proline"
                else:
                    # probably not a proline
                    print "not proline"

        else:
            # Not a residue so four clique carbon is either in a side chain or a terminal residue.
            # Why "hwel"?
            print "hwel"

        return [back_bone_list, main_chain_list, side_chain_list, residue_list]

    def draw_graph(self, G, node_list=None, edge_colour='k', node_size=15, node_colour='r', graph_type='spring',
                   back_bone=None, side_chains=None, terminators=None):
        # determine nodelist
        if node_list is None:
            node_list = G.nodes()
        # determine labels
        labels = {}
        for l_atom in G.nodes_iter():
            labels[l_atom] = l_atom.symbol

        # draw graphs based on graph_type
        if graph_type == 'circular':
            nx.draw_circular(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                             edge_color=edge_colour, node_color=node_colour)
        elif graph_type == 'random':
            nx.draw_random(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                           edge_color=edge_colour, node_color=node_colour)
        elif graph_type == 'spectral':
            nx.draw_spectral(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                             edge_color=edge_colour, node_color=node_colour)
        elif graph_type == 'spring':
            nx.draw_spring(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                           edge_color=edge_colour, node_color=node_colour)
        elif graph_type == 'shell':
            nx.draw_shell(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                          edge_color=edge_colour, node_color=node_colour)
        # elif graph_type == 'protein':
        # self.draw_protein(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
        #                   edge_color=edge_colour, node_color=node_colour, back_bone, side_chains, terminators)
        else:
            nx.draw_networkx(G, with_labels=True, labels=labels, node_list=node_list, node_size=node_size,
                             edge_color=edge_colour, node_color=node_colour)
        plt.show()

    # def draw_protein(self, G, back_bone, side_chains, terminators, num_turns, edge_colour='k', node_size=15,
    #                  node_colour='r'):
    # ''' This plotting routine assumes the graph is a protein.
    # It identifies the backbone and side chains.
    # It then plots the backbone on a spiral and computes sensible positions for the side chains.
    # '''
    # # set the length scale (essentially length of each edge in x - y coords)
    # length_scale = 1
    #
    # # find the longest side chains - sets a spacing parameter for the backbone spiral.
    # len_side_chains = [ nx.diameter(side_chain) for side_chain in side_chains]
    # max_side_chain = max(len_side_chains)
    #
    # # find the overall length (number of nodes in the back_bone)
    # len_back_bone = len(back_bone)
    #
    # # set angle per node: designed to give T full turns. 2 * Pi * T / (Num Nodes)
    # apn = 2.0 * np.pi * float(num_turns) / (float(len_back_bone - 1))
    #
    # # radial increase per node. Designed to give 2 * max_side_chain spacing when y = 0
    # rpn = 2.0 * float(num_turns) * float(max_side_chain) / (float(len_back_bone - 1))
    #
    # # compute the back bone spiral positions as a look up dictionary in back_bone_position
    # back_bone_pos = {}
    # for s, n in enumerate(back_bone):
    # back_bone_pos[n] = (rpn * s * length_scale * np.cos(apn * s), rpn * s * length_scale * np.sin(apn * s))
    #
    #
    # # identify the base positions of each side chain (c alpha positions in back bone)
    # for n in back_bone:
    # side_chain_base_position = {}
    #
    # nx.draw_networkx(G, with_labels=True, labels=labels, node_size=node_size, edge_color=edge_colour,
    #                  node_color=node_colour)

    def identify_chains(self):
        """ Creates sub graphs of the disconnected graphs in self.topology (i.e protein chains)
            and assigns a unique letter to them between A and Z.
            If there are more than 26 chains then the labeling beings at A again and so may not be unique."""

        # Create lists of the atoms in each disconnected sub component of the graph.
        chain_list = nx.connected_components(self.topology)

        # convert the list into a lists of sub graphs which "point" to the original graph.
        # nx.connected_component_subgraphs(self.topology) produces copies of the original graph.
        self.chains = []
        for c in chain_list:
            self.chains.append(self.topology.subgraph(c))

        # store the number of chains
        self.number_of_chains = len(self.chains)

        # Assign each of the sublists a unique letter.
        # On the day we find a protein with more than 26 chains we'll worry about that then.
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.chain_names = [alphabet[chain_id % 26] for chain_id in range(0, self.number_of_chains)]

    def find_longest_paths(self):
        """ finds the longest path in each chain. """

        self.longest_paths = []

        # loop through the chains
        for chain in self.chains:
            # The longest of the shortest paths between
            # the set of nodes in the periphery of the network is identified.

            # The eccentricity of a node is the maximum distance from a node to the other nodes.
            # the periphery is the set of nodes with the maximum eccentricity.
            # this is known as the diameter is network
            periphery = nx.periphery(chain)

            # find all the atoms belonging to the shortest paths between the elements of periphery 
            candidate_longest_paths = [
                [nx.shortest_path(chain, source, target) for target in periphery[source_index + 1:]] for
                source_index, source in enumerate(periphery)]

            # flatten each longest path
            candidate_longest_paths = [item for sub_list in candidate_longest_paths for item in sub_list]

            # get the length of each longest path
            longest_path_lengths = [len(long_path) for long_path in candidate_longest_paths]

            # find the max length (many will be the same, doesn't matter. get the first one).
            longest_path_index = np.argmax(longest_path_lengths)

            # extract the longest path from the current list of longest paths.
            # get one longest path per chain.
            self.longest_paths.append(candidate_longest_paths[longest_path_index])

    def find_cliques(self):
        """ Function finds the number of cliques to which each vertex belongs.
            Generates a separate list for each chain."""
        self.cliques = []
        for chain in self.chains:
            # for each chain generate a dictionary of the number of cliques 
            # to which each atom in the chain belongs
            self.cliques.append(nx.number_of_cliques(chain))

    def find_four_clique_carbons(self):
        """ Analyse each chain to determine the four clique carbons
            which also belong to the longest chain.
            Here's the logic:
            Carbons can belong to four cliques because they have four bonds.
            Typically all atoms on the longest path are C, N or H.
            N and H cannot belong to four cliques.
            In a protein the C=O of the peptide bond is always a double bond which means that is only has three cliques.
            Thus all four clique vertices on the longest path are either carbon alphas, carbons in the terminating
            groups, or carbons in the side chain of the terminating residue. But this list will contain ALL carbon
            alphas and, possibly, at most four extras per chain.
        """
        self.four_clique_carbons = []
        for chain, longest_path, clique in zip(self.chains, self.longest_paths, self.cliques):
            self.four_clique_carbons.append([atom for atom in chain if atom in longest_path and clique[atom] == 4])
