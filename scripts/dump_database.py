from optparse import OptionParser

def main():
    # add some program options
    parser = OptionParser(usage = "usage: %prog [options] storage")
    
    from pygmin.storage.database import Database
    
    parser.add_option("--write-disconnect",
                      dest="writeDPS", action="store_true",
                      help="generate min.dat and ts.dat to use with disconnectDPS")
    parser.add_option("-m",
                      dest="writeMinima", action="store_true",
                      help="dump minima to screen")
    parser.add_option("-t",
                      dest="writeTS", action="store_true",
                      help="dump transition states to screen")
    parser.add_option("-d",
                      dest="write_distances", action="store_true",
                      help="dump distances to screen")
    
    (options, args) = parser.parse_args()
    
    # print help if no input file is given
    if(len(args) != 1):
        parser.print_help()
        exit(-1)
        
    db = Database(db=args[0])

    if(options.writeMinima):
        print "List of minima:"
        print "---------------"
        for m in db.minima():
            print "%f\t\tid %d"%(m.energy, m._id)
        print "END\n"
    
    if(options.writeTS):
        print "List of transition states:"
        print "--------------------------"
        for ts in db.transition_states():
            print "%d\t<->\t%d\tid %d\tenergies %f %f %f"%\
                (ts.minimum1._id, ts.minimum2._id, ts._id, ts.minimum1.energy, ts.energy, ts.minimum2.energy)
        print "END\n"
    if(options.write_distances):
        print "List of distances:"
        print "--------------------------"
        for d in db.distances():
            print "%d\t<->\t%d\tid %d\tdistance %f"%\
                (d._minimum1_id, d._minimum2_id, d._id, d.dist)
        print "END\n"

    if(options.writeDPS):
        writeDPS(db)
        

def writeDPS(db):
    minindex={}
    out = open("min.data", "w")
    i=1
    for m in db.minima():
        minindex[m]=i
        i+=1
        if m.fvib is None:
            print "no frequency information available for minimum", m._id
        if m.pgorder is None:
            print "no pgorder information available for minimum", m._id
        out.write("%f %f %d 0.0 0.0 0.0\n"%(m.energy, m.fvib, m.pgorder))
    out = open("ts.data", "w")
    ti=0
    for ts in db.transition_states():
        ti+=1
        out.write("%f %f %d %d %d 0.0 0.0 0.0\n"%(ts.energy, ts.fvib, ts.pgorder, minindex[ts.minimum1], minindex[ts.minimum2]))
    print "Written %d minima and %d transition states"%(i, ti)

if __name__ == "__main__":
    main()    