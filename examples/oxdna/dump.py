import os
from optparse import OptionParser
import oxdnagmin_ as GMIN
from pygmin.storage.database import Database

TO_PDB="python /home/vr274/opt/oxDNA/UTILS/traj2vis.py  pdb %s gmindnatop"

def main():
    # add some program options
    parser = OptionParser(usage = "usage: %prog [options] storage")
        
    parser.add_option("--write-disconnect",
                      dest="writeDPS", action="store_true",
                      help="generate min.dat and ts.dat to use with disconnectDPS")
    parser.add_option("-m",
                      dest="writeMinima", action="store_true",
                      help="dump minima to screen")
    parser.add_option("-t",
                      dest="writeTS", action="store_true",
                      help="dump transition states to screen")
    parser.add_option("--coords",
                  dest="writeCoords", action="store_true",
                  help="export coordinates files")

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
        
    if(options.writeDPS):
        writeDPS(db)
        
    if(options.writeCoords):
        GMIN.initialize()
        i=0
        for m in db.minima():
            i+=1
            filename = "lowest/lowest%03d.cif"%(i)
            print "minimum",i, "energy",m.energy,"to",filename
            GMIN.userpot_dump(filename, m.coords)            
            if(not TO_PDB is None):
                os.system(TO_PDB%filename)

def writeDPS(db):
    minindex={}
    out = open("min.data", "w")
    i=1
    for m in db.minima():
        minindex[m]=i
        i+=1
        out.write("%f 0.0 1 0.0 0.0 0.0\n"%(m.energy))
    out = open("ts.data", "w")
    ti=0
    for ts in db.transition_states():
        ti+=1
        out.write("%f 0.0 1 %d %d 0.0 0.0 0.0\n"%(ts.energy, minindex[ts.minimum1], minindex[ts.minimum2]))
    print "Written %d minima and %d transition states"%(i, ti)

if __name__ == "__main__":
    main()    