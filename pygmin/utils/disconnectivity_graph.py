import networkx as nx
from collections import deque

from pygmin.landscape import Graph


        
class Tree(object):
    """
    a Tree graph
    """
    def __init__(self, parent=None):
        self.subtrees = []
        self.data = {}
        self.parent=parent
    
    def make_branch(self):
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
        if len(self.subtrees) == 0:
            nleaves = 1
        else:
            nleaves = 0
            for tree in self.subtrees:
                nleaves += tree.number_of_leaves()
        return nleaves
    
    def get_leaves(self):
        if len(self.subtrees) == 0:
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
        
        >>> graphwrapper = Graph(database)
        >>> dg = DisconnectivityGraph(graphwrapper.graph)
         
    nlevels : int
        how many levels at which to bin the transition states
    subgraph_size : int
        if subgraph_size is not None then all disconnected graphs
        of size greater than subraph_size will be included.
        
    """
    def __init__(self, graph, min0=None, nlevels=20,
                 subgraph_size=None):
        self.graph = graph
        self.nbins = nlevels
        self.subgraph_size = subgraph_size
        self.min0 = min0
        if self.min0 is None:
            #find the minimum energy node
            elist = [ (m.energy, m) for m in self.graph.nodes() ]
            elist = sorted(elist)
            self.min0 = elist[0][1]
            print "min0", self.min0.energy, self.min0._id
        self.transition_states = nx.get_edge_attributes(self.graph, "ts")
    
    def _getTS(self, min1, min2):
        try:
            return self.transition_states[(min1, min2)]
        except KeyError:
            return self.transition_states[(min2, min1)]
            
    
#    def _get_path(self, min1, min2):
#        #this should be a weighted path
#        return nx.shortest_path(self.graph, min1, min2)
#    
#    def _get_max_ts(self, min1, min2):
#        path = self._get_path(min1, min2)
#        emax = -1e100
#        for i in range(len(path)-1):
#            m1, m2 = path[i], path[i+1]
#            try:
#                edge = (m1, m2)
#                ts = self.transition_states[edge]
#            except KeyError:
#                edge = (m2, m1)
#                ts = self.transition_states[edge]
#
#            if ts.energy > emax:
#                tsmax = ts
#                emax = ts.energy
#
#        return tsmax

    def _connected_component_subgraphs(self, G):
        """
        redefine the networkx version because they use deepcopy
        on the nodes and edges
        
        this was coppied and only slightly modified from the original
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
            if self._getTS(*edge).energy < ethresh:
                newgraph.add_edge(*edge)
#        return nx.connected_component_subgraphs(newgraph)
        return self._connected_component_subgraphs(newgraph)


    def _recursive_new_tree(self, parent_tree, graph, energy_levels, ilevel):
        if ilevel < 0: return
        ethresh = energy_levels[ilevel]
#        if graph.number_of_edges() == 0: return
        
#        leaves = {}
        subgraphs = self._split_graph(graph, ethresh)
        for g in subgraphs:
            if nx.number_of_nodes(g) == 1:
                for minimum in g.nodes(): break
                newtree = parent_tree.make_branch()
                newtree.data["minimum"] = minimum 
                newtree.data["ilevel"] = ilevel 
                newtree.data["ethresh"] = ethresh 
#                leaves[minimum] = newtree
            else:
                newtree = parent_tree.make_branch()
                newtree.data["ilevel"] = ilevel 
                newtree.data["ethresh"] = ethresh                 
                self._recursive_new_tree(newtree, g, energy_levels, ilevel-1)

    def _make_tree(self, graph, energy_levels):
        """
        """
        tree_graph = Tree()
#        print tree_graph
        tree_graph.data["ilevel"] = len(energy_levels)-1
        de = energy_levels[-1] - energy_levels[-2]
        tree_graph.data["ethresh"] = energy_levels[-1] + 2*de
        
    
        #deal with the case that we have multiple disconnected graphs
        subgraphs = self._connected_component_subgraphs(graph)
        if len(subgraphs) > 1:
            #if the highest level tree has more than one graph then 
            #they are disconnected graphs
            tree_graph.data["children_not_connected"] = True
        
        for g in subgraphs:
            newtree = tree_graph.make_branch()
            newtree.data["ilevel"] = len(energy_levels)-1
            newtree.data["ethresh"] = energy_levels[-1] + de
            self._recursive_new_tree(newtree, g, energy_levels, len(energy_levels)-1)

        return tree_graph

#    def _recursive_new_tree_new(self, parent_tree, graphlist, energy_levels, ilevel):
#        """
#        for each graph in graphlist, make a branch on tree, split
#        graph according to energy_levels and pass on
#        """
#        if ilevel < 0: return
#        ethresh = energy_levels[ilevel]
#        
#        if len(graphlist) == 1:
#            #this is a leaf
#            g = graphlist[0]
#            for minimum in g.nodes(): break
#            newtree = parent_tree.make_branch()
#            newtree.data["minimum"] = minimum 
#            newtree.data["ilevel"] = ilevel 
#            newtree.data["ethresh"] = ethresh 
#
#        else:
#            for g in graphlist:
#                #split g
#                subgraphs = self._split_graph(g, ethresh)
#                newtree = parent_tree.make_branch()
#                newtree.data["ilevel"] = ilevel 
#                newtree.data["ethresh"] = ethresh                 
#                self._recursive_new_tree(newtree, subgraphs, energy_levels, ilevel-1)
                
            
        


    def _recursive_layout_x_axis(self, tree, xmin, dx_per_min):
#        nbranches = tree.number_of_branches()
        nminima = tree.number_of_leaves()
        subtrees = tree.get_subtrees()
        subtrees = self._order_trees_by_global_minimum(subtrees)
        tree.data["x"] = xmin + dx_per_min * nminima / 2.
        x = xmin
        for subtree in subtrees:
            self._recursive_layout_x_axis(subtree, x, dx_per_min)
            nminima_sub = subtree.number_of_leaves()
            x += dx_per_min * nminima_sub
  
    def _layout_x_axis(self, tree):
        """
        """
        xmin = 0.
        dx_per_min = 1.
        self._recursive_layout_x_axis(tree, xmin, dx_per_min)

    def _tree_get_minimum_energy(self, tree, emin=1e100):
        if tree.is_leaf():
            energy = tree.data["minimum"].energy
            if energy < emin:
                emin = energy
        else:
            for subtree in tree.get_subtrees():
                emin = self._tree_get_minimum_energy(subtree, emin)
        return emin

    def _order_trees_by_global_minimum(self, trees):
        """
        order trees with by the lowest energy minimum.  the global minimum
        goes in the center, with the remaining being placed alternating on the
        left and on the right.
        """
        mylist = [ (self._tree_get_minimum_energy(tree), tree) for tree in trees]
        mylist = sorted(mylist)
        neworder = deque()
        for i in range(len(mylist)):
            if i % 2 == 0:
                neworder.append(mylist[i][1])
            else:
                neworder.appendleft(mylist[i][1])
        return list(neworder)
            
        
# this is roughly the beginning of the algorithm in disconnectionDPS.f90    
#    def _assign_minima_to_basins(self, energy_levels):
#        """ energy_levels must be sorted from lowest to highest"""
#        
#        basin = {}
#        for n in self.graph.nodes(): basin[n] = -1
#        
#        count = 0
#        nbasin = 0
#        for ethresh in energy_levels:
#            changed  = True
#            while changed:
#                changed = False
#                if True:
#                    count += 1
#                    print "looping through again", count
#                for edge in self.graph.edges():
#                    ts = self.transition_states[edge]
#                    if ts.energy < ethresh:
#                        m1, m2 = ts.minimum1, ts.minimum2
#                        basin1, basin2 = basin[m1], basin[m2]
##                        print "%d\t%d\t%d\t%d" % (m1._id, m2._id, basin1, basin2)
#                        if basin1 == -1 and basin2 == -1:
#                            #if both minima are unasigned, assigne them to nbasin
#                            changed = True
#                            basin[m1] = nbasin
#                            basin[m2] = nbasin
#                            nbasin += 1
#                        elif basin1 != basin2:
#                            #they were previously assigned different basins
#                            changed = True
#                            if basin1 == -1:
#                                basin[m1] = basin2
#                            elif basin2 == -1:
#                                basin[m2] = basin1
#                            else:
#                                #assign them both to the same basin
#                                #use the lower number because that is the higher energy basin
#                                basin[m1] = min(basin1, basin2)
#                                basin[m1] = min(basin1, basin2)
#                        
#        return basin
                        

    def _get_energy_levels(self, graph):
        #define the energy levels
        elist = [self._getTS(*edge).energy for edge in graph.edges()]
        emin = min(elist)
        emax = max(elist)
        de = (emax - emin) / (self.nbins-1)
        #the lower edge of the bins
        elower = [emin + de*i for i in range(self.nbins)]
        
        return elower

    def _get_line_segment_recursive(self, line_segments, tree, eoffset):
        """
        add the line segment connecting this tree to it's parent
        """
        is_leaf = tree.number_of_leaves() == 1
        if tree.parent is not None:
            xparent = tree.parent.data['x']
            x = tree.data['x']
            yparent = tree.parent.data["ethresh"]
            if is_leaf:
                ylow = tree.data["minimum"].energy
            else:
                ylow = tree.data["ethresh"]
            yhigh = yparent - eoffset
            #add vertical line segment
            line_segments.append( ([x,x], [ylow, yhigh]) )
            if not tree.parent.data.has_key("children_not_connected"):
                #add angled line segment
                line_segments.append( ([xparent, x], [yparent,yhigh]) )
        for subtree in tree.get_subtrees():
            self._get_line_segment_recursive(line_segments, subtree, eoffset)

        
    def _get_line_segments(self, tree, eoffset=1.):
        """
        get all the line segments for drawing the connection between 
        each minimum to it's parent node.
        """
        line_segments = []
        self._get_line_segment_recursive(line_segments, tree, eoffset)
        return line_segments
    
    def get_minima_layout(self):
        """
        return the x position of the minima        
        """
        leaves = self.tree_graph.get_leaves()
        minima = [leaf.data["minimum"] for leaf in leaves]
        xpos = [leaf.data["x"] for leaf in leaves]
        return xpos, minima
    
    def _reduce_graph(self, graph, min0):
        """determine how much of the graph
        to include in the disconnectivity graph
        """
        used_nodes = []
        #make sure we include the subgraph containing min0
        used_nodes += nx.node_connected_component(self.graph, self.min0)
        
        if self.subgraph_size is not None:
            node_lists = nx.connected_components(self.graph)
            for nodes in node_lists:
                if len(nodes) >= self.subgraph_size:
                    used_nodes += nodes

        graph = self.graph.subgraph(used_nodes)
        return graph

        
    
    def calculate(self):
        """
        do the calculations necessary to draw the diconnectivity graph
        """
        #find a reduced graph with only those connected to min0
#        nodes = nx.node_connected_component(self.graph, self.min0)
#        self.graph = self.graph.subgraph(nodes)
        graph = self._reduce_graph(self.graph, self.min0)

        #define the energy levels
        elower = self._get_energy_levels(graph)
        
        #make the tree graph defining the discontinuity of the minima
        tree_graph = self._make_tree(graph, elower)
        
        #layout the x positions of the minima and the nodes
        self._layout_x_axis(tree_graph)

        #get the line segments which will be drawn to define the graph
        eoffset = (elower[-1] - elower[-2]) * 0.2  #this should be passable
        line_segments = self._get_line_segments(tree_graph, eoffset=eoffset)
        
        self.tree_graph = tree_graph
        self.line_segments = line_segments
    
    def plot(self, show_minima=False):
        #draw the minima as points
        import matplotlib.pyplot as plt
        leaves = self.tree_graph.get_leaves()
        energies = [leaf.data["minimum"].energy for leaf in leaves]
        xpos = [leaf.data["x"] for leaf in leaves]
        
        if show_minima:      
            plt.plot(xpos, energies, 'o')
        
        #draw the line segemnts
        for x, y in self.line_segments:
            plt.plot(x, y, 'k')
            
        plt.ylabel("Energy")
        plt.xticks([])
        plt.box(on=False)

        
    
        
#        import pylab as pl
#        pl.ioff()
##        nx.draw(self.graph, with_labels=False)
##        pl.show()
#        nx.draw(tree_graph, with_labels=False)
#        pl.show()
#        print "done drawing"
        
        
        
    