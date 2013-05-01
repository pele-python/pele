import networkx as nx
from collections import deque
from pygmin.landscape import Graph

__all__ = ["DisconnectivityGraph", "database2graph"]

def database2graph(database):
    """create a networkx graph from a pygmin database"""
    graph_wrapper = Graph(database)
    return graph_wrapper.graph
    

class Tree(object):
    """
    a Tree graph
    
    Each member of this class is a node in a Tree.  The node can
    have many children, but only one parent.  If the node has no
    parents then it is the root node.  If the node has no children
    then it is a leaf.
            
    """
    def __init__(self, parent=None):
        self.subtrees = []
        self.data = {}
        self.parent=parent
    
    def make_branch(self):
        """return a new Tree which is a child of this Tree"""
        newtree = Tree(parent=self)
        self.subtrees.append(newtree)
        return newtree
    
    def get_subtrees(self):
        return self.subtrees
    
    def number_of_branches(self):
        return len(self.subtrees)
    
    def is_leaf(self):
        return self.number_of_branches() == 0
    
    def number_of_leaves(self):
        """return the number of leaves that are descendants of this Tree""" 
        if len(self.subtrees) == 0:
            nleaves = 1
        else:
            nleaves = 0
            for tree in self.subtrees:
                nleaves += tree.number_of_leaves()
        return nleaves
    
    def get_leaves(self):
        """return a list of the leaves that are descendants of this Tree""" 
        if self.is_leaf():
            leaves = [self]
        else:
            leaves = []
            for tree in self.subtrees:
                leaves += tree.get_leaves()
        return leaves

                

class DisconnectivityGraph(object):
    """
    make a disconnectivity graph
    
    Parameters
    ----------
    graph : a networkx graph
        a graph with Minimum objects as nodes and transition
        states defining the edges.  You can use
        pygmin.landscape.Graph to create this from a database.
        
        >>> from pygmin.landscape import Graph
        >>> graphwrapper = Graph(database)
        >>> dg = DisconnectivityGraph(graphwrapper.graph)
         
    nlevels : int
        how many levels at which to bin the transition states
    Emax : float
        maximum energy for transition state cutoff.  Default
        is the maximum energy of all transition states.

    minima : list of Minima
        a list of minima to ensure are displayed in the graph.
        e.g. they will be displayed even if they're in a connected
        cluster smaller than smaller than subgraph_size
    subgraph_size : int
        if subgraph_size is not None then all disconnected graphs
        of size greater than subraph_size will be included.

    order_by_energy : bool
        order the subtrees by placing ones with lower energy closer
        to the center
    order_by_basin_size : bool
        order the subtrees by placing larger basins closer to the center
    center_gmin : bool
        when a node splits into its daughter
        nodes, the one containing the global minimum is always placed centrally
        (even if other nodes carry more minima). This does not guarantee that
        the global minimum is central in the overall diagram because other
        nodes may push the one containing the global minimum over to one side
    include_gmin : bool
        make sure to include the global minimum, even if it is not part of the
        main connected region
    node_offset : float
        offset between 0 and 1 for how to draw the angled lines.
        0 for no angle, draw horizontally out then vertically down.  1 for 
        angled lines all the way to the next energy level        
    energy_attribute : string, optional
        attribute which contains energy. default is energy. This attribute can
        be used to generate free energy disconnectivity graphs
    
    See Also
    ---------
    make_disconnectivity_graph.py :
        a script (in pygmin/scripts) to make the disconnectivity graph from the command line
    pygmin.storage.Database :
        The database format in which minima and transition states are stored in pygmin
    pygmin.landscape.Graph : 
        a wrapper to create a networkx Graph from a database
    
    Examples
    --------
    These examples assume a Database with minima already exists
    
    >>> from pygmin.landscape import Graph
    >>> import matplotlib.pyplot as plt
    >>> graphwrapper = Graph(database)
    >>> dg = DisconnectivityGraph(graphwrapper.graph)
    >>> dg.calculate()
    >>> dg.plot()
    >>> plt.show()

    """
    def __init__(self, graph, minima=None, nlevels=20, Emax=None,
                 subgraph_size=None, order_by_energy=False,
                 order_by_basin_size=True, node_offset=1.0,
                 center_gmin=True, include_gmin=True, energy_attribute="energy"):
        self.graph = graph
        self.nlevels = nlevels
        self.Emax = Emax
        self.subgraph_size = subgraph_size
        self.order_by_basin_size = order_by_basin_size
        self.order_by_energy = order_by_energy
        self.center_gmin = center_gmin 
        self.gmin0 = None
        self.energy_attribute = energy_attribute
        self.node_offset = node_offset
        if self.center_gmin:
            include_gmin = True

        if minima is None:
            minima = []
        self.min0list = minima
        if include_gmin:
            #find the minimum energy node
            elist = [ (self._getEnergy(m), m) for m in self.graph.nodes() ]
            elist = sorted(elist)
            self.gmin0 = elist[0][1]
            self.min0list.append(self.gmin0)
#            print "min0", self.min0.energy, self.min0._id
        self.transition_states = nx.get_edge_attributes(self.graph, "ts")
        self.minimum_to_leave = dict()
        self.tree_list = [[] for x in range(self.nlevels)]

    def _getEnergy(self, node):
        """ get the energy of a node """
        return getattr(node, self.energy_attribute)
    
    def set_energy_levels(self, elevels):
        ''' manually set the energy levels '''
        self.elevels = elevels
        
    def _getTS(self, min1, min2):
        """return the transition state object between two minima"""
        try:
            return self.transition_states[(min1, min2)]
        except KeyError:
            return self.transition_states[(min2, min1)]


    #############################################################
    #functions for building the tree by splitting the graph into
    #connected components at each level 
    #############################################################


    def _connected_component_subgraphs(self, G):
        """
        redefine the networkx version because they use deepcopy
        on the nodes and edges
        
        this was copied and only slightly modified from the original
        source
        """
        cc=nx.connected_components(G)
        graph_list=[]
        for c in cc:
            graph_list.append(G.subgraph(c))
        return graph_list
    
    def _split_graph(self, graph, ethresh):
        """
        make a new graph from the old one after excluding edges with energy > ethresh
        
        return a list of graphs, one for each connected component.
        """
        newgraph = nx.Graph()
        newgraph.add_nodes_from(graph.nodes())
        for edge in graph.edges():
            if self._getEnergy(self._getTS(*edge)) < ethresh:
                newgraph.add_edge(*edge)
#        return nx.connected_component_subgraphs(newgraph)
        return self._connected_component_subgraphs(newgraph)


    def _recursive_new_tree(self, parent_tree, graph, energy_levels, ilevel):
        """the recursive function called by `_make_tree` to build the disconnectivity Tree"""
        if ilevel < 0: return #  this the past the lowest level.  We're done
        ethresh = energy_levels[ilevel]
#        if graph.number_of_edges() == 0: return
        
#        leaves = {}
        subgraphs = self._split_graph(graph, ethresh)
        for g in subgraphs:
            if nx.number_of_nodes(g) == 1:
                #  There's only one mimimum in this graph.  Attach the minimum
                #  object to the tree node.
                for minimum in g.nodes(): break
                newtree = parent_tree.make_branch()
                newtree.data["minimum"] = minimum 
                newtree.data["ilevel"] = ilevel 
                newtree.data["ethresh"] = ethresh 
                self.minimum_to_leave[minimum] = newtree
#                leaves[minimum] = newtree
            else:
                #  subdivide this graph at the next energy level down.
                newtree = parent_tree.make_branch()
                newtree.data["ilevel"] = ilevel 
                newtree.data["ethresh"] = ethresh                 
                self._recursive_new_tree(newtree, g, energy_levels, ilevel-1)

    def _make_tree(self, graph, energy_levels):
        """make the disconnectivity graph tree
        
        start at the highest energy level, and at each energy level Elevel, remove the 
        edges with energy higher than Elevel.  This breaks the graph into disconnected
        components (subgraphs), which become nodes in the disconnectivity tree.
        Recursively repeat the process for each of those subgraphs in order to 
        build the disconnectivity graph.
        """
        tree_graph = Tree()
#        print tree_graph
        tree_graph.data["ilevel"] = len(energy_levels)-1
        de = energy_levels[-1] - energy_levels[-2]
        tree_graph.data["ethresh"] = energy_levels[-1] + 1.*de
        
    
        #deal with the case that we have multiple disconnected graphs
#        subgraphs = self._connected_component_subgraphs(graph)
        subgraphs = self._split_graph(graph, energy_levels[-1])
        if len(subgraphs) > 1:
            #if the highest level tree has more than one graph then 
            #they are disconnected graphs
            tree_graph.data["children_not_connected"] = True
        
        for g in subgraphs:
            newtree = tree_graph.make_branch()
            newtree.data["ilevel"] = len(energy_levels)-2
            newtree.data["ethresh"] = energy_levels[-1]
            self._recursive_new_tree(newtree, g, energy_levels, len(energy_levels)-2)

        return tree_graph            

    ##########################################################
    #functions for determining the x position of the branches
    #and leaves
    ##########################################################
    
    def _recursive_assign_id(self, tree):
        subtrees = tree.get_subtrees()
        for subtree in subtrees:
            if subtree.number_of_branches() >= 2:
                self.tree_list[subtree.data['ilevel']].append(subtree)

                subtree.data['id'] = len(self.tree_list[subtree.data['ilevel']])

            self._recursive_assign_id(subtree)

            
    def _assign_id(self, tree):
        """
        Determining the id of the branches and leaves
        for selection purposes
        """

        self._recursive_assign_id(tree)
        
    def _set_colour(self,i,colour_dict):
        '''
        
        '''
        self.tree_list[i[0]][i[1]].data['colour'] = colour_dict[i]

    def assign_colour(self, tree, colour_dict=[]):
        '''
        Colour trees according to `colour_dict`, a dictionay with 
        (level, tree_index) tuples as keys and RGB colours as values
        '''
        for i in colour_dict: self._set_colour(i,colour_dict)
            
            

    def _recursive_layout_x_axis(self, tree, xmin, dx_per_min):
#        nbranches = tree.number_of_branches()
        nminima = tree.number_of_leaves()
        subtrees = tree.get_subtrees()
        subtrees = self._order_trees(subtrees)
        tree.data["x"] = xmin + dx_per_min * nminima / 2.
        x = xmin
        for subtree in subtrees:
            self._recursive_layout_x_axis(subtree, x, dx_per_min)
            nminima_sub = subtree.number_of_leaves()
            x += dx_per_min * nminima_sub
  
    def _layout_x_axis(self, tree):
        """determining the x position of the branches and leaves
        
        used in displaying the disconnectivity graph
        """
        xmin = 4.0
        dx_per_min = 1
        self._recursive_layout_x_axis(tree, xmin, dx_per_min)

    def _tree_get_minimum_energy(self, tree, emin=1e100):
        """
        return the minimum energy of all the leaves in the tree
        """
        return min([leaf.data["minimum"].energy for leaf in tree.get_leaves()])

    def _order_trees(self, trees):
        """
        order a list of trees for printing
        
        This is the highest level function for ordering trees.  This, 
        and functions called by this, will account for all the user options 
        like center_gmin and order by energy
        """
        if self.order_by_energy:
            return self._order_trees_by_minimum_energy(trees)
        else:
            return self._order_trees_by_most_leaves(trees)

    def _order_trees_final(self, tree_value_list):
        """
        Parameters
        ----------
        tree_value_list :
            a list of (value, tree) pairs where value is the object
            by which to sort the trees
        
        Returns
        -------
        treelist :
            a list of trees ordered with the lowest in the center
            and the others placed successively on the left and right
        """
        mylist = sorted(tree_value_list)
        neworder = deque()
        for i in range(len(mylist)):
            if i % 2 == 0:
                neworder.append(mylist[i][1])
            else:
                neworder.appendleft(mylist[i][1])
        return list(neworder)

    def _ensure_gmin_is_center(self, tree_value_list):
        """ensure that the tree containing the global minimum has the lowest value
        """
        if self.gmin0 is None: return
        min0index = None
        for i in range(len(tree_value_list)):
            v, tree = tree_value_list[i]
            if self.gmin0 in [ leaf.data["minimum"] for leaf in tree.get_leaves()]:
                min0index = i
                break
        if min0index is not None:
            minvalue = min([v for v, tree in tree_value_list])
            #replace the value with a lower one
            #for the tree containing min0
            newvalue = minvalue - 1 #this won't work for non number values
            tree_value_list[i] = (newvalue, tree_value_list[i][1]) 
        return tree_value_list 

            
        

    def _order_trees_by_most_leaves(self, trees):
        """order list of trees by the number of leaves"""
#        if self.center_gmin:
#            return self._order_trees_by_most_leaves_and_global_min(trees)
        mylist = [ (tree.number_of_leaves(), tree) for tree in trees]
        if self.center_gmin:
            mylist = self._ensure_gmin_is_center(mylist)
        return self._order_trees_final(mylist)


    def _order_trees_by_minimum_energy(self, trees):
        """
        order trees with by the lowest energy minimum.  the global minimum
        goes in the center, with the remaining being placed alternating on the
        left and on the right.
        """
        mylist = [ (self._tree_get_minimum_energy(tree), tree) for tree in trees]
        return self._order_trees_final(mylist)
                        

    #######################################################################
    #functions which return the line segments that make up the visual graph
    #######################################################################

    def _get_line_segment_recursive(self, line_segments,line_colours, tree, eoffset):
        """
        add the line segment connecting this tree to it's parent
        """
        is_leaf = tree.is_leaf() == 1
        if tree.parent is not None:
            xparent = tree.parent.data['x']
            x = tree.data['x']
            yparent = tree.parent.data["ethresh"]
            
            try: line_colours.append(tree.parent.data['colour'])
            except KeyError: line_colours.append((0.0,0.0,0.0))
                     
            if is_leaf:
                ylow = self._getEnergy(tree.data["minimum"])
            else:
                ylow = tree.data["ethresh"]
            
            yhigh = max(yparent - eoffset, ylow)
                        
            if(yparent - eoffset > ylow or tree.number_of_branches() > 0):
                #add vertical line segment
                line_segments.append( ([x,x], [ylow, yhigh]) )
            else: # stop diagonal line earlier to avoid artifacts
                x = (x-xparent)/eoffset * (yparent - ylow) + xparent
                
            if not tree.parent.data.has_key("children_not_connected"):
                #add angled line segment
                line_segments.append( ([xparent, x], [yparent,yhigh]) )
                
        for subtree in tree.get_subtrees():
            self._get_line_segment_recursive(line_segments, line_colours, subtree, eoffset)

        
    def _get_line_segments(self, tree, eoffset=-1.):
        """
        get all the line segments for drawing the connection between 
        each minimum to it's parent node.
        """
        line_segments = []
        line_colours = []
        self._get_line_segment_recursive(line_segments, line_colours, tree, eoffset)
        return line_segments, line_colours
    
    
    ##########################################################################
    # functions for determining which minima to include in the 
    # disconnectivity graph
    ##########################################################################
    
    def _remove_nodes_with_few_edges(self, graph, nmin):
        rmlist = [n for n in graph.nodes() if graph.degree(n) < nmin]
        if len(rmlist) > 0:
            if self.gmin0 is not None:
                if self.gmin0 in rmlist:
                    print "global minimum has", graph.degree(self.gmin0), "edges, not showing in graph"
            print "removing", len(rmlist), "minima from graph with fewer than", nmin, "edges"
            for n in rmlist:
                graph.remove_node(n)
        return graph

    
    def _remove_high_energy_minima(self, graph, emax):
        if emax is None: return graph
        rmlist = [m for m in graph.nodes() if self._getEnergy(m) > emax]
        if len(rmlist) > 0:
            print "removing %d nodes with energy higher than"%len(rmlist), emax
        for m in rmlist:
            graph.remove_node(m)
        return graph

    def _remove_high_energy_transitions(self, graph, emax):
        if emax is None: return graph
        rmlist = [edge for edge in graph.edges() \
                  if self._getEnergy(self._getTS(edge[0], edge[1])) > emax]
        if len(rmlist) > 0:
            print "removing %d edges with energy higher than"%len(rmlist), emax
        for edge in rmlist:
            graph.remove_edge(edge[0], edge[1])
        return graph

    def _reduce_graph(self, graph, min0list):
        """determine how much of the graph to include in the disconnectivity graph
        """
        used_nodes = []
        #make sure we include the subgraph containing min0
        if len(min0list) == 0:
            #use the biggest connected cluster
            cc = nx.connected_components(graph)
            used_nodes += cc[0] #list is ordered by size of cluster
        else:
            for min0 in min0list:
                used_nodes += nx.node_connected_component(graph, min0)
        
        if self.subgraph_size is not None:
            node_lists = nx.connected_components(graph)
            for nodes in node_lists:
                if len(nodes) >= self.subgraph_size:
                    used_nodes += nodes

        newgraph = graph.subgraph(used_nodes)
        return newgraph

    ##########################################################################
    # general functions
    ##########################################################################
    
    
    def get_minima_layout(self):
        """
        return the x position of the minima        
        """
        leaves = self.tree_graph.get_leaves()
        minima = [leaf.data["minimum"] for leaf in leaves]
        xpos = [leaf.data["x"] for leaf in leaves]
        return xpos, minima
    
    def get_tree_layout(self):
        '''
        Returns the x position of the trees
        ''' 
        id = []

        for l in range(len(self.tree_list)):
            id += [tuple([l,i]) for i in range(len(self.tree_list[l]))]
        x_pos = [self.tree_list[l][i].data['x'] for l, i in id] 
        energies = [self.tree_list[l][i].data['ethresh'] for l, i in id]

        return id, x_pos, energies
        
    def _get_energy_levels(self, graph):
        """
        combine input and the graph data to determine what the 
        energy levels will be.
        """
        
        if hasattr(self, "elevels"):
            return self.elevels
        
        #define the energy levels
        elist = [self._getEnergy(self._getTS(*edge)) for edge in graph.edges()]
        if len(elist) == 0:
            raise Exception("there are no edges in the graph.  Is the global minimum connected?")
        emin = min(elist)
        if self.Emax is None:
            emax = max(elist)
        else:
            emax = self.Emax
        de = (emax - emin) / (self.nlevels-1)
        #the upper edge of the bins
        elower = [emin + de*(i) for i in range(self.nlevels)]
        elevels = elower.append(emin + de*self.nlevels)
        
        return elower

       
    
    def calculate(self):
        """
        do the calculations necessary to draw the diconnectivity graph
        """
        graph = self.graph
        assert graph.number_of_nodes() > 0, "graph has no minima"
        assert graph.number_of_edges() > 0, "graph has no transition states"
        
        # we start with applying the energy cutoff, otherwise reduce
        # graph does not work as intended
        graph = self._remove_high_energy_minima(graph, self.Emax)
        graph = self._remove_high_energy_transitions(graph, self.Emax)
        assert graph.number_of_nodes() > 0, "after applying Emax, graph has no minima"
        assert graph.number_of_edges() > 0, "after applying Emax, graph has no minima" 
        
        #find a reduced graph with only those connected to min0
#        nodes = nx.node_connected_component(self.graph, self.min0)
#        self.graph = self.graph.subgraph(nodes)
        graph = self._reduce_graph(graph, self.min0list)
        
        #define the energy levels
        elevels = self._get_energy_levels(graph)
        
        #remove more nodes
        graph = self._remove_high_energy_minima(graph, elevels[-1])
        graph = self._remove_high_energy_transitions(graph, elevels[-1])
        graph = self._remove_nodes_with_few_edges(graph, 1)
        
        assert graph.number_of_nodes() > 0, "after cleaning up the graph, graph has no minima"
        assert graph.number_of_edges() > 0, "after cleaning up the graph, graph has no minima" 

        #make the tree graph defining the discontinuity of the minima
        tree_graph = self._make_tree(graph, elevels)
        
        #assign id to trees
        self._assign_id(tree_graph)
        
        #assign colour to trees
        self.assign_colour(tree_graph)#, colour_list)
        
        #layout the x positions of the minima and the nodes
        self._layout_x_axis(tree_graph)

        #get the line segments which will be drawn to define the graph
        eoffset = (elevels[-1] - elevels[-2]) * self.node_offset  #this should be passable
        line_segments = self._get_line_segments(tree_graph, eoffset=eoffset)
        
        self.eoffset = eoffset
        self.tree_graph = tree_graph
        self.line_segments = line_segments
    
    def plot(self, show_minima=False, show_trees=False, linewidth=0.5, axes=None):
        """draw the disconnectivity graph using matplotlib
        
        don't forget to call calculate() first
        
        also, you must call pyplot.show() to actually see the plot
        """
        import matplotlib as mpl
        from matplotlib.collections import LineCollection
        import matplotlib.pyplot as plt
        
        self.line_segments, self.line_colours = self._get_line_segments(self.tree_graph, eoffset=self.eoffset)
        
        #set up how the figure should look
        if axes is not None:
            ax = axes
        else:
            fig = plt.figure(figsize=(6,7))
            fig.set_facecolor('white')
            ax = fig.add_subplot(111, adjustable='box')

        ax.tick_params(axis='y', direction='out')
        ax.yaxis.tick_left()
        ax.spines['left'].set_color('black')
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['right'].set_color('none')
#        plt.box(on=True)

#         if show_trees:
#             trees = self.tree_graph.get_subtrees()
        
        #draw the minima as points
        if show_minima:      
            leaves = self.tree_graph.get_leaves()
            energies = [self._getEnergy(leaf.data["minimum"]) for leaf in leaves]
            xpos = [leaf.data["x"] for leaf in leaves]
        
            ax.plot(xpos, energies, 'o')
        
        # draw the line segments 
        # use LineCollection because it's much faster than drawing the lines individually 
        linecollection = LineCollection([ [(x[0],y[0]), (x[1],y[1])] for x,y in self.line_segments])
        linecollection.set_linewidth(linewidth)
        linecollection.set_color(self.line_colours)
        ax.add_collection(linecollection)
        
        # scale the axes appropriately
        ax.relim()
        ax.autoscale_view(scalex=True, scaley=True, tight=None)

        #remove xtics            
        ax.set_xticks([])        
        
        
    
