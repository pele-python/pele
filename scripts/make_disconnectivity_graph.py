import sys
import os
import getopt
import time
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import pele.utils.disconnectivity_graph as dg
from pele.storage import Database
from pele.utils.optim_compatibility import OptimDBConverter


try:
    from PyQt4.QtGui import QApplication
    from pele.gui.ui.dgraph_dlg import DGraphDialog, reduced_db2graph
    use_gui = True
except ImportError:
    use_gui = False

def read_minA(fname, db):
    """load data from min.A or min.B"""
    with open(fname) as fin:
        ids = []
        for i, line in enumerate(fin):
            if i == 0:
                nminima = int(line.split()[0])
            else:
                sline = line.split()
                ids += map(int, sline)
    
    assert nminima == len(ids)
    print len(ids), "minima read from file:", fname
    return [db.getMinimum(mid) for mid in ids]

def read_AB(db):
    minA = "min.A"
    minB = "min.B"
    if os.path.isfile(minA) and os.path.isfile(minA):
        A = read_minA(minA, db) 
        B = read_minA(minB, db)
        return A,B
    else:
        return None

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
    
    
    groups = None
    
    if OPTIM:
        #make database from min.data ts.data
        db = Database()
        converter = OptimDBConverter(db)
        converter.convert_no_coords()
        groups = read_AB(db)
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
        if groups is not None:
            md.dgraph_widget.dg.color_by_group(groups)
            md.dgraph_widget.redraw_disconnectivity_graph()
        md.show()
        sys.exit(app.exec_())
        
    # make graph from database
    t0 = time.time()
    if "Emax" in kwargs and use_gui:
        graph = reduced_db2graph(db, kwargs['Emax'])
    else:
        graph = dg.database2graph(db)
    t1 = time.time()
    print "loading the data into a transition state graph took", t1-t0, "seconds"

    # do the disconnectivity graph analysis
    mydg = dg.DisconnectivityGraph(graph, **kwargs)
    print "doing disconnectivity graph analysis"
    sys.stdout.flush()
    t1 = time.time()
    mydg.calculate()
    t2 = time.time()
    print "d-graph analysis finished in", t2-t1, "seconds"
    print "number of minima:", mydg.tree_graph.number_of_leaves()
    print "plotting disconnectivigy graph"
    sys.stdout.flush()
    
    
    # make the figure and save it
    mydg.plot()
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)
    t3 = time.time()
    print "plotting finished in", t3-t2, "seconds"
        
    

if __name__ == "__main__":
    main()