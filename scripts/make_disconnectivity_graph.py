import sys
import os
import getopt
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import pygmin.utils.disconnectivity_graph as dg
from pygmin.storage import Database


try:
    from PyQt4.QtGui import QApplication
    from pygmin.gui.ui.dgraph_dlg import DGraphDialog, reduced_db2graph
    use_gui = True
except ImportError:
    use_gui = False
    

global _id_count
_id_count = 0

class FakeMinimum(object):
    """
    a class to duplicate some of the functionality of the Minimum class
    """
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
        global _id_count
        self._id = _id_count
        _id_count += 1

    def __eq__(self, m):
        """m can be integer or Minima object"""
        assert self._id is not None
        if isinstance(m, FakeMinimum):
            assert m._id is not None
            return self._id == m._id
        else:
            return self._id == m
        
    def __hash__(self):
        assert self._id is not None
        return self._id
    

class FakeTransitionState(object):
    """
    a class to duplicate some of the functionality of the 
    TransitionState class
    """
    def __init__(self, energy, coords, min1, min2):
        self.energy = energy
        self.coords = coords
        self.minimum1 = min1
        self.minimum2 = min2
        global _id_count
        self._id = _id_count
        _id_count += 1


def make_graph_from_files(mindata="min.data", tsdata="ts.data"):
    """skip making the database because it takes too long"""
    graph = nx.Graph()
    id2min = {}
    print "adding minima from", mindata
    count = 0
    with open(mindata, "r") as fin:
        for line in fin:
            sline = line.split()
            energy = float(sline[0])
            coords = np.array([0.])
            minimum = FakeMinimum(energy, coords)
            id2min[minimum._id] = minimum
            graph.add_node(minimum)
#            count += 1
#            if count % 1000 == 0:
#                print count, minimum._id

    count = 0
    print "adding transition states from", tsdata
    with open(tsdata, "r") as fin:
        for line in fin:
            sline = line.split()
            energy = float(sline[0])
            m1id = int(sline[3]) - 1
            m2id = int(sline[4]) - 1
            m1 = id2min[m1id]
            m2 = id2min[m2id]
            coords = np.array([0.])
            ts = FakeTransitionState(energy, coords, m1, m2)
            graph.add_edge(m1, m2, ts=ts)
#            count += 1
#            if count % 100 == 0:
#                print count

    return graph 


def usage():
    print "usage:"
    print sys.argv[0], "database [options]"
    print "   database is the file which contains your database"
    print "   -o outfile : save the plot to this pdf file."
    print "   --OPTIM :    load data from min.data and ts.data"
    print ""
    print " options to pass to DisconnectivityGraph:"
    print "   --nlevels=n : number of energy levels"
    print "   --subgraph_size=n :  include all disconnected subgraphs up to size n"
    print "   --order_by_basin_size : " 
    print "   --order_by_energy : " 
    print "   --include_gmin : " 
    print "   --center_gmin : "
    print "   --Emax=emax   :" 

def main():
    if len(sys.argv) < 2:
        usage()
        exit(1)
    
    
    
    kwargs = {}
    outfile = None
    OPTIM = False

    opts, args = getopt.gnu_getopt(sys.argv[1:], "ho:", 
                                   ["help", "nlevels=", "subgraph_size=", "OPTIM",
                                    "order_by_basin_size", "order_by_energy",
                                    "include_gmin",
                                    "center_gmin",
                                    "Emax=",
                                    ])
    for o, a in opts:
        if o == "-h" or o == "--help":
            usage()
            exit(1)
        if o == "-o":
            outfile = a 
        elif o == "--nlevels":
            kwargs["nlevels"] = int(a)
        elif o == "--Emax":
            kwargs["Emax"] = float(a)
        elif o == "--subgraph_size":
            kwargs["subgraph_size"] = int(a)
        elif o == "--order_by_basin_size":
            kwargs["order_by_basin_size"] = True
        elif o == "--order_by_energy":
            kwargs["order_by_energy"] = True
        elif o == "--includer_gmin":
            kwargs["includer_gmin"] = True
        elif o == "--center_gmin":
            kwargs["center_gmin"] = True
        elif o == "--OPTIM":
            OPTIM = True
        else:
            print "don't understand", o, a
            print ""
            usage()
            exit(1)
    
    
    if OPTIM:
        #make database from min.data ts.data
        graph = make_graph_from_files()
    else:
        if len(args) == 0:
            print "you must specify database file"
            print ""
            usage()
            exit()
        dbfile = args[0]
        if not os.path.exists(dbfile):
            print "database file doesn't exist", dbfile
            exit()
        
        db = Database(dbfile)
        
    if outfile is None and use_gui:
        app = QApplication(sys.argv) 
        kwargs["show_minima"] = False
        md = DGraphDialog(db, params=kwargs)
        md.rebuild_disconnectivity_graph()
        md.show()
        sys.exit(app.exec_())
        
    if not OPTIM: graph = reduced_db2graph(db, kwargs['Emax'])
    mydg = dg.DisconnectivityGraph(graph, **kwargs)
    mydg.calculate()
    print "number of minima:", mydg.tree_graph.number_of_leaves()
    mydg.plot()
    
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)
        
    

if __name__ == "__main__":
    main()