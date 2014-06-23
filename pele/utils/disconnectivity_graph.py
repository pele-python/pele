import copy
from itertools import izip
from collections import deque

import numpy as np
import networkx as nx

from pele.landscape import database2graph

__all__ = ["DisconnectivityGraph"]

class TreeLeastCommonAncestor(object):
    """Find the least common ancestor to a set of trees"""
    def __init__(self, trees):
        self.start_trees = trees
        self.run()
    
    def run(self):
        # find all common ancestors
        common_ancestors = set()
        for tree in self.start_trees:
            parents = set(tree.get_ancestors())
            parents.add(tree)
            if len(common_ancestors) == 0:
                common_ancestors.update(parents)
            else:
                # remove all elements that are not common
                common_ancestors.intersection_update(parents)
                assert len(common_ancestors) > 0

        if len(common_ancestors) == 0:
            raise Exception("the trees don't have any common ancestors")
        
        # sort the common ancestors by the number of ancestors each has
        common_ancestors = list(common_ancestors)
        if len(common_ancestors) > 1:
            common_ancestors.sort(key=lambda tree: len(list(tree.get_ancestors())))

        # the least common ancestor is the one with the most ancestors
        self.least_common_ancestor = common_ancestors[-1]
        return self.least_common_ancestor
    
    def get_all_paths_to_common_ancestor(self):
        """return all the ancestors of all the input trees up to the least common ancestor"""
        trees = set(self.start_trees)
        for tree in self.start_trees:
            for parent in tree.get_ancestors():
                trees.add(parent)
                if parent == self.least_common_ancestor:
                    break
        return trees
        
            
#            for tree in common_ancestors:
#                for parent in tree.get_ancestors():
#                    if parent in common_ancestors
#            
#        return iter(common_ancestors).next()


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
        self.parent = parent
    
    def add_branch(self, branch):
        """make branch a child of this tree""" 
        self.subtrees.append(branch)
        branch.parent = self
    
    def make_branch(self):
        """return a new Tree which is a child of this Tree"""
        newtree = self.__class__(parent=self)
        self.subtrees.append(newtree)
        return newtree
    
    def get_subtrees(self):
        print "get_subtrees() is deprecated. use get_branches() instead!"
        return self.get_branches()
    
    def get_branches(self):
        """return the list of branches of this tree"""
        return self.subtrees
    
    def number_of_branches(self):
        """return the number of branches of this tree"""
        return len(self.subtrees)
    
    def is_leaf(self):
        """return true if this tree has no descendants"""
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
    
    def leaf_iterator(self):
        """iterate through the leaves that are descendants of this Tree""" 
        if self.is_leaf():
            yield self
        else:
            for tree in self.subtrees:
                for leaf in tree.leaf_iterator():
                    yield leaf
    
    def get_all_trees(self):
        """iterator over all subtrees, including self"""
        yield self
        for branch in self.get_branches():
            for subtree in branch.get_all_trees():
                yield subtree

    def number_of_subtrees(self):
        ntot = 1
        for branch in self.get_branches():
            ntot += branch.number_of_subtrees()
        return ntot
    
    def get_ancestors(self):
        """iterate over ancestors excluding self"""
        if self.parent is not None:
            yield self.parent
            for ancestor in self.parent.get_ancestors():
                yield ancestor
            

class DGTree(Tree):
    """add a few functions to Tree to make it specific to disconnectivity graph"""
    def contains_minimum(self, min1):
        for leaf in self.get_leaves():
            if leaf.data["minimum"] == min1:
                return True
        return False
    
    def get_minima(self):
        return [leaf.data["minimum"] for leaf in self.get_leaves()]
    
    def get_one_minimum(self):
        """return a single minimum that is in this tree"""
        key = "_random_minimum"
        if self.is_leaf():
            return self.data["minimum"]
        elif key in self.data:
            return self.data[key]
        else:
            m = self.get_branches()[0].get_one_minimum()
            self.data[key] = m
            return m
    
    def _test_tree(self):
        tset = set()
        for tree in self.get_all_trees():
            if tree in tset:
                print "tree is touched twice"
                return False
            tset.add(tree)
        return True
        
class UnionFind(nx.utils.UnionFind):
    def groups_iter(self):
        return (c for c, c1 in self.parents.iteritems() if c == c1)

class _MakeTree(object):
    """class to Make the disconnectivity graph tree
    
    Parameters
    ----------
    minima : list of Minimum objects
        the minima to contain in the disconnectivity graph
    transition_states : list of TransitionState objects
        the transition states to contain in the disconnectivity graph
    energy_levels : list of floats
        the energy levels at which to split the graph
    get_energy : callable
        a function which returns the energy of a transition state.  You can 
        use this to redefine the `energy` of a transition state, e.g. for a 
        free energy disconnectivity graph.
        
    Notes
    -----
    This algorithm starts from a completely disconnected graph and adds
    transition states one at a time, starting from the lowest energy transition
    state.  Connectivity is determined solely through labeling minima based on
    which cluster it's in.  If a new transition
    state connects two previously disconnected clusters, the clusters are
    joined to make one single cluster.  As the transition states are added
    (sorted in energy) the energy levels are reached one at a time.  At each
    level, the state of the connectivity of the graph is saved in tree graphs.
    
    This algorithm is very similar to kruskal's minimum spanning tree algorithm
    
    """
    def __init__(self, minima, transition_states, energy_levels, get_energy=None):
        self.minima = minima
        self.transition_states = transition_states
        self.energy_levels = energy_levels
        self._get_energy = get_energy
        
        self.union_find = UnionFind()
        self.minimum_to_leave = dict()


    def get_energy(self, ts):
        """return the energy of the transition state"""
        if self._get_energy is None:
            return ts.energy
        else:
            return self._get_energy(ts)

    def _new_leaf(self, m):
        """make a new leaf from minimum m"""
        leaf = DGTree()
        leaf.data["minimum"] = m 
        self.minimum_to_leave[m] = leaf
        return leaf
        
    def make_tree(self):
        """make the disconnectivity tree"""
        # make list of transition states sorted so that lower energies are to the right
        tslist = filter(lambda ts: ts.minimum1 != ts.minimum2, self.transition_states)
        # remove duplicate entries and sort
        tslist = list(set(tslist))
        tslist.sort(key=lambda ts: -self.get_energy(ts))
        self.transition_states = tslist 

        trees = []
        
        # build the tree up starting at the lowest level
        for ilevel in range(len(self.energy_levels)):
            trees = self._do_next_level(ilevel, trees)

        # deal with any disconnected parts
        energy_levels = self.energy_levels
        if len(trees) == 1:
            self.tree = trees[0]
        else:
            self.tree = DGTree()
            for t in trees:
                self.tree.add_branch(t)
            self.tree.data["ilevel"] = len(energy_levels)-1
            de = energy_levels[-1] - energy_levels[-2]
            self.tree.data["ethresh"] = energy_levels[-1] + 1.*de
            self.tree.data["children_not_connected"] = True

        if False:
            res = self.tree._test_tree()
            print "tree test result", res

        return self.tree

    def _add_edge(self, min1, min2):
        """add an edge between min1 and min2
        
        if min1 and min2 belong to different color groups, those groups
        will be set equal
        """
        new_minima = []
        if min1 not in self.union_find.parents:
            new_minima.append(min1)
        if min2 not in self.union_find.parents:
            new_minima.append(min2)
        self.union_find.union(min1, min2)
        return new_minima

    def _do_next_level(self, ilevel, previous_trees):
        """do the disconnectivity analysis for energy level ilevel
        """
        ethresh = self.energy_levels[ilevel]
        
        # add the edges to the graph up to ethresh
        tslist = self.transition_states
        while len(tslist) > 0:
            ts = tslist[-1]
            if self.get_energy(ts) >= ethresh:
                break
            
            new_minima = self._add_edge(ts.minimum1, ts.minimum2)
            for m in new_minima:
                previous_trees.append(self._new_leaf(m))
            
            # remove the transition state from the list
            tslist.pop()
        
        # make a new tree for every color (connected cluster)
        newtrees = []
        color_to_tree = dict()
        for c in self.union_find.groups_iter():
            newtree = DGTree()
            newtree.data["ilevel"] = ilevel 
            newtree.data["ethresh"] = ethresh
            newtrees.append(newtree)
            color_to_tree[c] = newtree
        
        # determine parentage
        for tree in previous_trees:
            m = tree.get_one_minimum()
            c = self.union_find[m]
            parent = color_to_tree[c]
            if tree.number_of_branches() == 1:
                # remove linear parentage.
                subtree = iter(tree.subtrees).next()
                parent.add_branch(subtree)
            else:
                parent.add_branch(tree)
            
        return newtrees


class ColorDGraphByGroups(object):
    """color the graph based on specified grouping of minima

    Parameters
    ----------
    tree_graph: a DGTree object
        usually accessed by dgraph.tree_graph if dgraph is a 
        DisconnectivityGraph object.
    groups : list
        list of groups of minima that should have the same color
    
    Notes
    -----
    For each node, check all minima for which the node is a parent.
    If all minima are contained on one of the groups, the node 
    will be coloured to represent that group.  
    If any minimum is not contained on one of the groups, the node
    is not coloured. 
    If all minima are contained in groups but more than one group 
    is represented, the node will be the colour of the last group listed
    """
    def __init__(self, tree_graph, groups, colors=None):
        self.tree_graph = tree_graph
        
        # set the colors
        self._minimum_to_color = dict()
        self.color_list = self.get_list_of_colors(len(groups), colors=colors)
        for color, group in izip(self.color_list, groups):
            for minimum in group:
                self._minimum_to_color[minimum] = color
        self._tree_to_colors = dict()
    
    def get_list_of_colors_mpl(self, number):
        """return a list of colors for the groups.  Use matplotlib colormap"""
        from matplotlib import cm
        colormap = cm.get_cmap("Dark2", lut=number)
        colors = [colormap(i) for i in np.linspace(0., 1., number)]
        return colors
    
    def parse_list_of_colors(self, number, colors):
        from matplotlib.colors import ColorConverter
        cconvert = ColorConverter()
        if number != len(colors):
            raise ValueError("the length of colors must be the number of groups")
        rgbcolors = [cconvert.to_rgb(c) for c in colors]
        return rgbcolors 
        
         
    def get_list_of_colors(self, number, colors=None):
        """return a list of colors for the groups"""
        if colors is not None:
            return self.parse_list_of_colors(number, colors)
        try:
            import brewer2mpl
        except ImportError:
            print "could not import brewer2mpl."
            print "install package brewer2mpl for a nicer color scheme"
            return self.get_list_of_colors_mpl(number)
        if number <= 12:
            bnumber = max(3, number)
            bcolors = brewer2mpl.get_map("Dark2", "Qualitative", bnumber)
            colors = bcolors.mpl_colors
            colors = colors[:number]
        else:
            return self.get_list_of_colors_mpl(number)
        return colors
    
    def minimum_to_color(self, minimum):
        """return the color of the minimum, or None if not colored"""
        try:
            return self._minimum_to_color[minimum]
        except KeyError:
            return None
    
    def tree_get_colors(self, tree):
        """return the color that this tree should be colored by"""
        try:
            return self._tree_to_colors[tree]
        except KeyError:
            if tree.is_leaf():
                color = self.minimum_to_color(tree.data["minimum"])
                if color is None:
                    colors = None
                else:
                    colors = frozenset([color])    
                self._tree_to_colors[tree] = colors
                return colors
            else:
                colors_list = [self.tree_get_colors(subtree) for subtree in tree.get_branches()]
                if None in colors_list:
                    colors = None
                else:
                    colors = frozenset([g for colors1 in colors_list for g in colors1])
                self._tree_to_colors[tree] = colors
                return colors
    
    def colors_to_color(self, colors):
        """if there are multiple colors for this tree, select which one to use
        
        select the color which is listed last in self.color_list
        """
        for color in reversed(self.color_list):
            if color in colors:
                return color       
            
    
    def run(self):
        """main loop for the algorithm"""
        for tree in self.tree_graph.get_all_trees():
            colors = self.tree_get_colors(tree)
            if colors is not None:
                tree.data["colour"] = self.colors_to_color(colors)
    
class ColorDGraphByValue(object):
    """color a disconnectivity graph by values associated with minima (e.g. order parameter)

    Parameters
    ----------
    tree_graph: a DGTree object
        usually accessed by dgraph.tree_graph if dgraph is a 
        DisconnectivityGraph object.
    minimum_to_value: callable
        A function that accepts a minimum and returns a float value.
        return None to indicate no color for this minimum
    colormap: callable, optional
        function which converts a float in (0,1) to a matplotlib color (RGB)
    normalize_values: bool
        if True the values will be normalized to fall between 0 and 1
    
    Notes
    -----
    Each node in the graph will be colored according to the value of the 
    child minimum with the largest value.  If any child minimum has value None
    then the node will not be colored
    """
    def __init__(self, tree_graph, minimum_to_value, colormap=None, 
                 normalize_values=True):
        self.tree_graph = tree_graph
        self.minimum_to_value = minimum_to_value
        if colormap is None:
            from matplotlib import cm
            self.colormap = cm.get_cmap("winter")
        else:
            self.colormap = colormap
    
        self._tree_to_value = dict()
        
        if normalize_values:
            values = [self.minimum_to_value(leaf.data["minimum"]) 
                      for leaf in self.tree_graph.leaf_iterator()]
            values = filter(lambda v: v is not None, values)
            self.maxval = max(values)
            self.minval = min(values)
        else:
            self.minval = None
            self.maxval = None
    
    def resolve_multiple_values(self, values):
        if None in values:
            return None
        else:
            return max(values)
    
    def value_to_color(self, value):
        if self.minval is None:
            vnorm = value
        else:
            vnorm = (value - self.minval) / (self.maxval - self.minval)
        return self.colormap(vnorm)
    
    def tree_get_value(self, tree):
        """return the color that this tree should be colored by"""
        try:
            return self._tree_to_value[tree]
        except KeyError:
            if tree.is_leaf():
                value = self.minimum_to_value(tree.data["minimum"])
                self._tree_to_value[tree] = value
                return value
            else:
                values = [self.tree_get_value(subtree) for subtree in tree.get_branches()]
                value = self.resolve_multiple_values(values)
                self._tree_to_value[tree] = value
                return value
    
            
    
    def run(self):
        """main loop for the algorithm"""
        for tree in self.tree_graph.get_all_trees():
            value = self.tree_get_value(tree)
            if value is not None:
                tree.data["colour"] = self.value_to_color(value)    
            
    
    

class DisconnectivityGraph(object):
    """
    make a disconnectivity graph
    
    Parameters
    ----------
    graph : a networkx graph
        a graph with Minimum objects as nodes and transition
        states defining the edges.  You can use the function
        database2graph() defined in this module to create this 
        from a database.
        
        >>> from pele.utils.disconnectivity_graph import database2graph
        >>> graph = database2graph(database)
        >>> dg = DisconnectivityGraph(graph)
         
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
    order_by_value : callable, optional, `v = order_by_value(m)`
        If this function is passed the minima at each level will be sorted by this value
        with small values to the left.  A group of minima will be sorted according to the 
        smallest value in the group.  
    
    See Also
    ---------
    make_disconnectivity_graph.py :
        a script (in pele/scripts) to make the disconnectivity graph from the command line
    pele.storage.Database :
        The database format in which minima and transition states are stored in pele
    
    Examples
    --------
    These examples assume a Database with minima already exists
    
    >>> import matplotlib.pyplot as plt
    >>> graph = database2graph(database)
    >>> dg = DisconnectivityGraph(graph)
    >>> dg.calculate()
    >>> dg.plot()
    >>> plt.show()

    """
    def __init__(self, graph, minima=None, nlevels=20, Emax=None,
                 subgraph_size=None, order_by_energy=False,
                 order_by_basin_size=True, node_offset=1.,
                 center_gmin=True, include_gmin=True, energy_attribute="energy",
                 order_by_value=None):
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
        self.get_value = order_by_value
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
    # functions for building the tree by splitting the graph into
    # connected components at each level 
    #############################################################


    def _connected_component_subgraphs(self, G):
        """
        redefine the networkx version because they use deepcopy
        on the nodes and edges
        
        this was copied and only slightly modified from the original
        source
        """
        cc = nx.connected_components(G)
        graph_list=[]
        for c in cc:
            graph_list.append(G.subgraph(c))
        return graph_list
    

    def _make_tree(self, graph, energy_levels):
        """make the disconnectivity graph tree
        """
        transition_states = nx.get_edge_attributes(graph, "ts").values()
        minima = graph.nodes()
        maketree = _MakeTree(minima, transition_states, energy_levels, 
                             get_energy=self._getEnergy)
        trees = maketree.make_tree()
        self.minimum_to_leave = maketree.minimum_to_leave
        return trees
        
    #################################################################
    # These functions determine how to layout the tree on the x axis
    #################################################################

    def _recursive_layout_x_axis(self, tree, xmin, dx_per_min):
#        nbranches = tree.number_of_branches()
        nminima = tree.number_of_leaves()
        subtrees = tree.get_branches()
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
        dx_per_min = 1.
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
        if self.get_value is not None:
            return self._order_trees_by_value(trees)
        elif self.order_by_energy:
            return self._order_trees_by_minimum_energy(trees)
        else:
            return self._order_trees_by_most_leaves(trees)

    def _order_trees_by_value(self, trees):
        """order the trees by a value. smaller numbers to the left
        
        Each tree will take the smallest value of all its associated minima.
        """
        def get_min_val(tree):
            return min([self.get_value(leaf.data["minimum"])
                        for leaf in tree.leaf_iterator()])
        trees.sort(key=get_min_val)
        return trees

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
    # functions which return the line segments that make up the visual graph
    #######################################################################

    def _get_line_segment_single(self, line_segments, line_colours, tree, eoffset):
        color_default = (0., 0., 0.)
        if tree.parent is None:
            # this is a top level tree.  Add a short decorative vertical line
            if "children_not_connected" in tree.data:
                # this tree is simply a container, add the line to the subtrees
                treelist = tree.subtrees
            else:
                treelist = [tree]
            
            dy = self.energy_levels[-1] - self.energy_levels[-2]
            for t in treelist:
                x = t.data["x"]
                y = t.data["ethresh"]
                line_segments.append(([x,x], [y,y+dy]))
                line_colours.append(color_default)
        else:
            # add two line segments.  A vertical one to yhigh
            #  ([x, x], [y, yhigh])
            # and an angled one connecting yhigh with the parent
            #  ([x, xparent], [yhigh, yparent])
            xparent = tree.parent.data['x']
            xself = tree.data['x']
            yparent = tree.parent.data["ethresh"]

            if tree.is_leaf():
                yself = self._getEnergy(tree.data["minimum"])
            else:
                yself = tree.data["ethresh"]
            
            # determine yhigh from eoffset
            yhigh = yparent - eoffset
            if yhigh <= yself:
                draw_vertical = False
                yhigh = yself
            else:
                draw_vertical = True
            
            # determine the line color
            try: 
                color = tree.data['colour']
            except KeyError:
                color = color_default
            
            # draw vertical line
            if tree.is_leaf() and not draw_vertical:
                # stop diagonal line earlier to avoid artifacts
                # change the x position so that the angle of the line
                # doesn't change
                if not "_x_updated" in tree.data:
                    dxdy = (xself - xparent) / eoffset
                    xself = dxdy * (yparent - yself) + xparent
                    tree.data['x'] = xself
                    tree.data["_x_updated"] = True
            else: 
                #add vertical line segment
                line_segments.append( ([xself,xself], [yself, yhigh]) )
                line_colours.append(color)
#                print "coloring vertical line", tree
            
            # draw the diagonal line
            if not tree.parent.data.has_key("children_not_connected"):
                line_segments.append( ([xself, xparent], [yhigh, yparent]) )
                line_colours.append(color)

    def _get_line_segment_recursive(self, line_segments, line_colours, tree, eoffset):
        """
        add the line segment connecting this tree to its parent
        """
        self._get_line_segment_single(line_segments, line_colours, tree, eoffset)
        for subtree in tree.get_branches():
            self._get_line_segment_recursive(line_segments, line_colours, subtree, eoffset)

        
    def _get_line_segments(self, tree, eoffset=-1.):
        """
        get all the line segments for drawing the connection between 
        each minimum to it's parent node.
        
        Returns
        -------
        line_segments : list
            list of line segments.  each line segment has the form ((x1, x2), (y1, y2))
        """
        line_segments = []
        line_colours = []
        self._get_line_segment_recursive(line_segments, line_colours, tree, eoffset)
        assert len(line_segments) == len(line_colours)
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
            cc = list(nx.connected_components(graph))
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
        
        # find a reduced graph with only those connected to min0
#        nodes = nx.node_connected_component(self.graph, self.min0)
#        self.graph = self.graph.subgraph(nodes)
        graph = self._reduce_graph(graph, self.min0list)
        
        # define the energy levels
        elevels = self._get_energy_levels(graph)
        self.energy_levels = elevels
        
        # remove more nodes
        graph = self._remove_high_energy_minima(graph, elevels[-1])
        graph = self._remove_high_energy_transitions(graph, elevels[-1])
        graph = self._remove_nodes_with_few_edges(graph, 1)
        
        assert graph.number_of_nodes() > 0, "after cleaning up the graph, graph has no minima"
        assert graph.number_of_edges() > 0, "after cleaning up the graph, graph has no edges" 

        # make the tree graph defining the discontinuity of the minima
        tree_graph = self._make_tree(graph, elevels)
        
        #assign id to trees
        # this is needed for coloring basins
#        self._assign_id(tree_graph)
        
        #layout the x positions of the minima and the nodes
        self._layout_x_axis(tree_graph)

        #get the line segments which will be drawn to define the graph
        eoffset = (elevels[-1] - elevels[-2]) * self.node_offset  #this should be passable
#        line_segments = self._get_line_segments(tree_graph, eoffset=eoffset)
        
        
        self.eoffset = eoffset
        self.tree_graph = tree_graph
#        self.line_segments = line_segments


    def color_by_group(self, groups, colors=None):
        """color the graph based on specified grouping of minima
    
        Parameters
        ----------
        groups : list
            list of groups of minima that should have the same color
        
        Notes
        -----
        For each node, check all minima for which the node is a parent.
        If all minima are contained on one of the groups, the node 
        will be coloured to represent that group.  
        If any minimum is not contained on one of the groups, the node
        is not coloured. 
        If all minima are contained in groups but more than one group 
        is represented, the node will be the colour of the last group listed
        """
        colorer = ColorDGraphByGroups(self.tree_graph, groups, colors=None)
        colorer.run()

    def color_by_value(self, minimum_to_value, colormap=None, 
                 normalize_values=True):        
        """color the graph by values associated with minima (e.g. order parameter)
    
        Parameters
        ----------
        minimum_to_value: callable
            A function that accepts a minimum and returns a float value.
            return None to indicate no color for this minimum
        colormap: callable, optional
            function which converts a float in (0,1) to a matplotlib color (RGB)
        normalize_values: bool
            if True the values will be normalized to fall between 0 and 1
        
        Notes
        -----
        Each node in the graph will be colored according to the value of the 
        child minimum with the largest value.  If any child minimum has value None
        then the node will not be colored
        """
        colorer = ColorDGraphByValue(self.tree_graph, minimum_to_value,
                                     colormap=colormap, 
                                     normalize_values=normalize_values)
        colorer.run()

    def plot(self, show_minima=False, show_trees=False, linewidth=0.5, axes=None,
             title=None):
        """draw the disconnectivity graph using matplotlib
        
        don't forget to call calculate() first
        
        also, you must call pyplot.show() to actually see the plot
        """
        from matplotlib.collections import LineCollection
        import matplotlib.pyplot as plt
        
        self.line_segments, self.line_colours = self._get_line_segments(self.tree_graph, eoffset=self.eoffset)

        # get the axes object        
        if axes is not None:
            ax = axes
        else:
            try:
                ax = self.axes
            except AttributeError:
                fig = plt.figure(figsize=(6,7))
                fig.set_facecolor('white')
                ax = fig.add_subplot(111, adjustable='box')

        #set up how the figure should look
        ax.tick_params(axis='y', direction='out')
        ax.yaxis.tick_left()
        # make the borders a bit prettier
        ax.spines['left'].set_color('black')
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['right'].set_color('none')
#        plt.box(on=True)

        if title is not None:
            ax.set_title(title)

        #draw the minima as points
        if show_minima: 
            xpos, minima = self.get_minima_layout()
            energies = [m.energy for m in minima]
            ax.plot(xpos, energies, 'o')
        
        # draw the line segments 
        # use LineCollection because it's much faster than drawing the lines individually 
        linecollection = LineCollection([ [(x[0],y[0]), (x[1],y[1])] for x,y in self.line_segments])
        linecollection.set_linewidth(linewidth)
        linecollection.set_color(self.line_colours)
        ax.add_collection(linecollection)
        
        
        # scale the axes appropriately
        # note: do not call ax.relim().  As of matplotlib version 1.3
        # ax.relim() does not take Collections into account so it does
        # not compute the limits correctly.  ax.autoscale_view() seems to work just fine
        ax.autoscale_view(scalex=True, scaley=True, tight=False)
        ax.set_ybound(upper=self.Emax)
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmin-0.5, xmax+0.5)
        
        # remove xtics
        # note: the xticks are removed after ax.autoscale_view() is called.
        # If it is the other way around the lines are too close the image border          
        ax.set_xticks([])
        self.axes = ax
        
    def label_minima(self, minima_labels, axes=None, 
                     rotation=60., **kwargs):
        """label the specified minima
        
        Parameters
        ----------
        minima_labels: dict
            dictionary with minima as keys and labels as values.
            i.e. label = minima_labels[minimum]
        axes: matplotlib axis object, optional
            The axes we are working on
        rotation: float
            angle (in degrees) of how much to rotate the text
        kwargs: kwargs
            additional keyword arguments are passed on to matplotlib 
            ax.set_xticklabels()
        
        Notes
        -----
        if the labels are outside of the figure bounding box you can fix it with
        plt.tight_layout() or fig.tight_layout() 
        """
        if axes is not None:
            ax = axes
        else:
            try:
                ax = self.axes
            except AttributeError:
                print "you must call plot() before label_minima()"
                raise
        leaves = filter(lambda leaf: leaf.data["minimum"] in minima_labels, 
                        self.tree_graph.leaf_iterator())
        xpos = [leaf.data["x"] for leaf in leaves]
        labels = [minima_labels[leaf.data["minimum"]] for leaf in leaves]
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=rotation, **kwargs)

#         rescale = True
#         if rescale:
#             import matplotlib.pyplot as plt
#             plt.tight_layout()
 
    
    def show(self):
        """simple wrapper for matplotlib.pyplot.show()"""
        from matplotlib import pyplot
        pyplot.show()
    def savefig(self, *args, **kwargs):
        """simple wrapper for matplotlib.pyplot.savefig()"""
        from matplotlib import pyplot
        pyplot.savefig(*args, **kwargs)
    
