import networkx as nx
from collections import deque

from pygmin.landscape import Graph


        
class Tree(object):
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
    """
    def __init__(self, graph, min0=None, nbins=20):
        self.graph = graph
        self.nbins = nbins
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
            
    
    def _get_path(self, min1, min2):
        #this should be a weighted path
        return nx.shortest_path(self.graph, min1, min2)
    
    def _get_max_ts(self, min1, min2):
        path = self._get_path(min1, min2)
        emax = -1e100
        for i in range(len(path)-1):
            m1, m2 = path[i], path[i+1]
            try:
                edge = (m1, m2)
                ts = self.transition_states[edge]
            except KeyError:
                edge = (m2, m1)
                ts = self.transition_states[edge]

            if ts.energy > emax:
                tsmax = ts
                emax = ts.energy

        return tsmax

       
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
        return nx.connected_component_subgraphs(newgraph)


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
        tree_graph.data["ethresh"] = energy_levels[-1] + de
        self._recursive_new_tree(tree_graph, graph, energy_levels, len(energy_levels)-1)

        return tree_graph
  
    def _recursive_layout_x_axis(self, tree, xmin, dx_per_min):
#        nbranches = tree.number_of_branches()
        nminima = tree.number_of_leaves()
        subtrees = tree.get_subtrees()
        tree.data["x"] = xmin + dx_per_min * nminima / 2.
        x = xmin
        for subtree in subtrees:
            self._recursive_layout_x_axis(subtree, x, dx_per_min)
            nminima_sub = subtree.number_of_leaves()
            x += dx_per_min * nminima_sub
            
        
        
          
    def _layout_x_axis(self, graph, tree):
        """
        """
        xmin = 0.
        dx_per_min = 1.
        self._recursive_layout_x_axis(tree, xmin, dx_per_min)
        
        
# this is roughly the algorithm in disconnectionDPS.f90    
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
                        

    def _get_energy_levels(self):
        
        elist = [self.transition_states[edge].energy for edge in self.graph.edges()]
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
            line_segments.append( ([x,x], [ylow, yhigh]) )
            line_segments.append( ([xparent, x], [yparent,yhigh]) )
        for subtree in tree.get_subtrees():
            self._get_line_segment_recursive(line_segments, subtree, eoffset)

        
    def _get_line_segments(self, tree, eoffset=1.):
        line_segments = []
        self._get_line_segment_recursive(line_segments, tree, eoffset)
        return line_segments
    
    def run(self):
        #find a reduced graph with only those connected to min0
        nodes = nx.node_connected_component(self.graph, self.min0)
        self.graph = self.graph.subgraph(nodes)

        #define the energy levels
        elower = self._get_energy_levels()
        
        #make the tree graph defining the discontinuity of the minima
        tree_graph = self._make_tree(self.graph, elower)
        
        print "minima", nx.number_of_nodes(self.graph)
        print "nodes", tree_graph.number_of_leaves()
#        print "edges", nx.number_of_edges(tree_graph)
        
        #layout the x positions of the minima and the nodes
        self._layout_x_axis(self.graph, tree_graph)

        #get the line segments which will be drawn to define the graph
        line_segments = self._get_line_segments(tree_graph, eoffset=1.)
        
        #draw the minima as points
        import matplotlib.pyplot as plt
        leaves = tree_graph.get_leaves()
        energies = [leaf.data["minimum"].energy for leaf in leaves]
        xpos = [leaf.data["x"] for leaf in leaves]        
        plt.plot(xpos, energies, 'o')
        
        #draw the line segemnts
        for x, y in line_segments:
            plt.plot(x, y, 'k')

        
        plt.show()
    
        
#        import pylab as pl
#        pl.ioff()
##        nx.draw(self.graph, with_labels=False)
##        pl.show()
#        nx.draw(tree_graph, with_labels=False)
#        pl.show()
#        print "done drawing"
        
        
        
    