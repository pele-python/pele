import argparse

import networkx as nx

from pele.storage.database import Database
from pele.utils.disconnectivity_graph import database2graph
from pele.utils.optim_compatibility import WritePathsampleDB

def print_system_properties(db, supress_long=True):
    if len(db.properties()) == 0: return
    print "System properties:"
    print "------------------"
    for p in db.properties():
        name, value = p.name(), p.value()
        str_value = str(value)
        if len(str_value) > 100 and supress_long:
            str_value = str_value[:80] + " '... output suppressed'"
        print "%10s:\t\t%s" % (name, str_value)
    print ""
        

def long_summary(db):
    nts = db.number_of_transition_states()
    if nts == 0:
        print "long summary not applicable: no transition states"
        return
    graph = database2graph(db)
    
    cclist = nx.connected_components(graph)
    
#    print "number of connected components", len(cclist)
    counts = dict()
    minimum_energy = dict()
    for cc in cclist:
        nc = len(cc)
        Emin = min((m.energy for m in cc))
        try:
            counts[nc] += 1
            if Emin < minimum_energy[nc]:
                minimum_energy[nc] = Emin 
        except KeyError:
            counts[nc] = 1
            minimum_energy[nc] = Emin
    counts = counts.items()
    counts.sort(key=lambda x:-x[0])
    print "Connectivity of the database:"
    for n, count in counts:
        if n == 1:
            print "%7d unconnected minima                : minimum energy = %s" % (count, minimum_energy[n])
        else:
            print "%7d connected clusters of size %7d: minimum energy = %s" % (count, n, minimum_energy[n]) 

def write_pathsample_db(db):
    writer = WritePathsampleDB(db)
    writer.write_db()

def main():
    parser = argparse.ArgumentParser(description="print information about the database")

    parser.add_argument("database", type=str, help="Database file name")
    
    parser.add_argument("--write-pathsample-db",
                      dest="write_pathsample", action="store_true",
                      help="generate a pathsample database by writing files min.data, ts.data, points.min, and points.ts")
    parser.add_argument("-m",
                      dest="writeMinima", action="store_true",
                      help="dump minima to screen")
    parser.add_argument("-t",
                      dest="writeTS", action="store_true",
                      help="dump transition states to screen")
    parser.add_argument("-p",
                      dest="properties", action="store_true",
                      help="print system properties")
    parser.add_argument("-s",
                      dest="summary", action="store_true",
                      help="print summary")
    parser.add_argument("-S",
                      dest="summary_long", action="store_true",
                      help="print long summary")
    args = parser.parse_args()
    
    if args.summary_long:
        args.summary = True
        
    db = Database(db=args.database, createdb=False)

    if args.properties or args.summary:
        print_system_properties(db)

    if args.summary:
        print "number of minima:", db.number_of_minima()
        print "number of transition states:", db.number_of_transition_states()
  
    if args.summary_long:
        long_summary(db)
        
    if args.writeMinima:
        print "List of minima: energy id fvib pgorder"
        print "---------------"
        for m in db.minima():
            print "%f\t\tid %d %s %s" % (m.energy, m._id, str(m.fvib), str(m.pgorder))
        print "END\n"
    
    if args.writeTS:
        print "List of transition states:"
        print "--------------------------"
        for ts in db.transition_states():
            print "%d\t<->\t%d\tid %d\tenergies %f %f %f" % \
                (ts.minimum1._id, ts.minimum2._id, ts._id, ts.minimum1.energy, ts.energy, ts.minimum2.energy)
        print "END\n"

    if args.write_pathsample:
        write_pathsample_db(db)
        


if __name__ == "__main__":
    main()    
