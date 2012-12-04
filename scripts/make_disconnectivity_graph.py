import sys
import os
import getopt
import matplotlib.pyplot as plt

import pygmin.utils.disconnectivity_graph as dg
from pygmin.storage import Database
from pygmin.landscape import Graph

def usage():
    print "usage:"
    print sys.argv[0], "database [-h --help -o outfile] [options]"
    print "   database is the file which contains your database"
    print "   options is a list of keyword pairs to be passed to"
    print "           DisconnectivityGraph using --keyword=value format"
    print "           e.g. --nlevels=20"
    print "   -o outfile : save the plot to this pdf file."

def main():
    if len(sys.argv) < 2:
        usage()
        exit(1)
    
    
    
    kwargs = {}
    outfile = None

    opts, args = getopt.gnu_getopt(sys.argv[1:], "ho:", 
                                   ["help", "nlevels=", "nts_min=", "subgraph_size="])
    for o, a in opts:
        if o == "-h" or o == "--help":
            usage()
            exit(1)
        if o == "-o":
            outfile = a 
        elif o == "--nlevels":
            kwargs["nlevels"] = int(a)
        elif o == "--nts_min":
            kwargs["nts_min"] = int(a)
        elif o == "--subgraph_size":
            kwargs["subgraph_size"] = int(a)
        else:
            print "don't understand", o, a
            usage()
            exit(1)
    
    
    dbfile = args[0]
    if not os.path.exists(dbfile):
        print "database file doesn't exists", dbfile
    
    
    
    db = Database(dbfile)
    graphwrapper = Graph(db)
    
    mydg = dg.DisconnectivityGraph(graphwrapper.graph, **kwargs)
    mydg.calculate()
    print "number of minima:", mydg.tree_graph.number_of_leaves()
    mydg.plot()
    
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)
    

if __name__ == "__main__":
    main()