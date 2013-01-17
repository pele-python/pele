import os
import subprocess
import numpy as np

from pygmin.optimize import Result


class PathInfoReader(object):
    """read path.info files"""
    def __init__(self, natoms, fname="path.info"):
        self.natoms = natoms
        self.fname = fname
    
    def get_next_line(self):
        return self.fin.readline()[:-1]
    
    def read(self):
        with open(self.fname, "r") as fin:
            self.fin = fin
            while True:
                try:
                    min1 = self.read_minimum(fin)
                    ts = self.read_minimum(fin)
                    min2 = self.read_minimum(fin)
                
                    yield min1, ts, min2
                except IndexError:
                    print "I think I'm done... ??"
                    break
            
    
    def read_coords(self, fin):
        coords = np.zeros([self.natoms, 3])
        for i in range(self.natoms):
            line = self.get_next_line()
            coords[i,:] = tuple(line.split()[:3])
        return coords.reshape(-1)
    
    def read_minimum(self, fin):
        """
        the minima and transition states are stored as
        
        1 line with energy
        1 line with point group info
        natoms lines with eigenvalues
        natoms lines with coords
        """
        res = Result()
        #read energy
        line = self.get_next_line()
#        print line
        res.energy = float(line.split()[0])

        #ignore the line with the point group
        line = self.get_next_line()
        
        res.eigenvalues = self.read_coords(fin)
        res.coords = self.read_coords(fin)
        
        return res
    
#    def read_ts(self, fin):
#        coords = self.read_coords(fin)
#        #read eigenvalues and eigenvectors?
#         
#        return coords

class SpawnOPTIM(object):
    """
    this class will control spawning of optim jobs and importing the results
    
    1. make temp directory
    #. make imput files (odata, finish, perm.allow, etc)
    #. call OPTIM
    #. when OPTIM finishes, load the results into the database
    """
    def __init__(self, coords1, coords2, OPTIM="OPTIM"):
        self.coords1 = coords1.reshape(-1,3)
        self.coords2 = coords2.reshape(-1,3)
        self.OPTIM = OPTIM
        
    def run(self):
        rundir = self.make_temporary_dir()
        self.rundir = rundir
        print "rundir", rundir
        self.make_input_files(rundir)
        self.call_optim(rundir, self.OPTIM)
        
    def make_input_files(self, rundir):
        
        #make odata file
        odata = rundir + "/odata"
        with open(odata, "w") as fout:
            self.write_odata(fout)
            self.write_odata_coords(self.coords1, fout)
        
        #make finish file
        finish = rundir + "/finish"
        with open(finish, "w") as fout:
            for xyz in self.coords2:
                fout.write( "%f %f %f\n" % tuple(xyz))

    
    def write_odata(self, fout):
        """write the odata file to the open file fout
        
        stop after the line POINTS.  i.e. don't write coordinates
        """
        raise NotImplementedError
    
    def write_odata_coords(self, coords, fout):
        """write the coords to odata
        
        need to label the atoms correctly for the system
        """
        raise NotImplementedError

    
    def make_temporary_dir(self):
        dname = "odata_spawn"
        if not os.path.isdir(dname):
            os.mkdir(dname)
        return dname
    
    def call_optim(self, rundir, OPTIM):
        curdir = os.getcwd()

        os.chdir(rundir)
        try:
            with open("Oout", "w") as fout:
                subprocess.call(OPTIM, stdout=fout)
        except OSError:
            print "exception raised in calling OPTIM:", OPTIM
            print "is the path correct?"
            os.chdir(curdir)
            raise
        os.chdir(curdir)
    
    def load_results(self, database):
        reader = PathInfoReader(self.sys.natoms, fname=self.rundir+"/path.info")
        newminima = set()
        newts = set()
        for min1res, tsres, min2res in reader.read():
#            print min1
            min1 = database.addMinimum(min1res.energy, min1res.coords)
            min2 = database.addMinimum(min2res.energy, min2res.coords)
            ts = database.addTransitionState(tsres.energy, tsres.coords, min1, min2)
            #I should probably get the eigenvector here
            print "adding ", min1._id, min2._id
#            print min1.energy, ts.energy, min2.energy
            newminima.add(min1)
            newminima.add(min2)
            newts.add(ts)
        
        return newminima, newts
        
    
class SpawnOPTIM_LJ(SpawnOPTIM):
    def __init__(self, coords1, coords2, sys, **kwargs):
        super(SpawnOPTIM_LJ, self).__init__(coords1, coords2, **kwargs)
        self.sys = sys
    
    def write_odata_coords(self, coords, fout):
        coords = coords.reshape(-1,3)
        for xyz in coords:
            fout.write( "AX %f %f %f\n" % tuple(xyz))
    
    def write_odata(self, fout):
        odatastr = """
DUMPALLPATHS

NEWCONNECT 50 1 15.0 14.0 22 5.0 0.025
ADJUSTK 5 10.0D0 1.05D0
USEDIAG 2
RANROT 10
REOPTIMISEENDPOINTS
NEBK 5.0
comment NEWNEB 50 500 0.0001
PERMDIST
DIJKSTRA 2
MAXBFGS  0.2 0.2
EDIFFTOL 0.000001
GEOMDIFFTOL 0.05
MAXSTEP 0.15
MAXMAX  0.125
TRAD 0.12
BFGSTS  50 3 25 0.0001
UPDATES 10 10
comment NOIT
NOHESS
BFGSMIN  0.000001
PUSHOFF 0.02
STEPS    200
BFGSSTEPS 5000
POINTS
"""
        fout.write(odatastr)
        fout.write("\n")

if __name__ == "__main__":
    from pygmin.systems import LJCluster
    natoms = 13
    sys = LJCluster(natoms)
    db = sys.create_database()
    x1, E1 = sys.get_random_minimized_configuration()[:2]
    x2, E2 = sys.get_random_minimized_configuration()[:2]
    m1 = db.addMinimum(E1, x1)
    m2 = db.addMinimum(E2, x2)
    
    optim = "/home/js850/git/OPTIM/source/build/OPTIM"
    spawner = SpawnOPTIM_LJ(x1, x2, sys, OPTIM=optim)
    spawner.run()
    spawner.load_results(db)
