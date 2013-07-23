import os
from optparse import OptionParser
import oxdnagmin_ as GMIN
from pele.storage.database import Database
import numpy as np
from pele.utils import rotations
from pele.utils.rbtools import CoordsAdapter
    
TO_PDB="python /home/vr274/opt/oxDNA/UTILS/traj2vis.py  pdb %s gmindnatop"

def export_xyz(fl, coords):
    ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
    fl.write("%d\n\n"%(2*ca.nrigid))
    for i in xrange(ca.nrigid):
        a = np.dot(rotations.aa2mx(ca.rotRigid[i]), np.array([1., 0., 0.]))
        x_back = ca.posRigid[i] - 0.4*a # backbone bead
        x_stack = ca.posRigid[i] + 0.4*a
        
        fl.write("C %f %f %f\n"%(x_back[0], x_back[1], x_back[2]))
        fl.write("H %f %f %f\n"%(x_stack[0], x_stack[1], x_stack[2]))

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
    parser.add_option("--xyz",
                  dest="writeXYZ", action="store_true",
                  help="export xyz files")

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
            np.savetxt("lowest/coords_%03d.txt"%(i), m.coords)
            
    if(options.writeXYZ):
        traj=open("lowest/traj.xyz", "w")
        i=0
        for m in db.minima():
            i+=1
            filename = "lowest/lowest%03d.xyz"%(i)
            print "minimum",i, "energy",m.energy,"to",filename
            export_xyz(open(filename, "w"), m.coords)
            export_xyz(traj, m.coords)

        traj.close()


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