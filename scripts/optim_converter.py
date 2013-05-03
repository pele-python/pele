import structure_read
import numpy as np
from pygmin.storage import Database, Minimum, TransitionState
# import sys
# import getopt
import argparse

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class Convert(object):
    '''
    Converts old OPTIM to pygmin database
    '''
    def __init__(self,db_name=None):
        if not db_name: db_name = "Database.db"
        self.db_name = db_name
        self.db = None
        self.natoms = None
        self.mindata = None
        self.tsdata = None
        self.pointsmin = None
        self.pointsts = None
                
    def setNatoms(self,natoms):
        self.natoms = natoms
        
    def setMindata(self,mindata=None):
        if not mindata: mindata = 'min.data'
        self.mindata = mindata
        
    def setTsdata(self,tsdata=None):
        if not tsdata: tsdata = 'ts.data'
        self.tsdata = tsdata
        
    def setPointsmin(self,pointsmin=None):
        if not pointsmin: pointsmin = 'points.min'
        self.pointsmin = pointsmin
        
    def setPointsts(self,pointsts=None):
        if not pointsts: pointsts = 'points.ts'
        self.pointsts = pointsts
        
    def initDatabase(self):
        self.db = Database(self.db_name)
        
    def setAccuracy(self,accuracy = 0.000001):
        self.db.accuracy = accuracy
        
    def ReadMindata(self):
        indx = 1
        f_len = file_len(self.mindata)
        for line in open(self.mindata,'r'):
            sline = line.split()
            coords = structure_read.structure_read(self.pointsmin,
                                                   indx,
                                                   self.natoms)
#             coords = np.random.random(3*self.natoms)
            e, f = map(float,sline[:2])
            pg = int(sline[2])
            
            min1 = Minimum(e, coords)
            
            min1.fvib = f
            min1.pgorder = pg
            
            min1.op = 1.0
            
            self.db.session.add(min1)

            indx += 1

    def ReadTSdata(self):
        indx = 1
        f_len = file_len(self.tsdata)
        for line in open(self.tsdata,'r'):
            sline = line.split()

            coords = structure_read.structure_read(self.pointsts,
                                                   indx,          
                                                   self.natoms)
#             coords = np.random.random(3*self.natoms)
            e, f = map(float,sline[:2])
            pg = int(sline[2])
            m1, m2 = map(int,sline[3:5])
            
            min1 = self.db.getMinimum(m1)
            min2 = self.db.getMinimum(m2)
            trans = TransitionState(e, coords, min1, min2)
            
            
            trans.fvib = f
            trans.pgorder = pg
            
            self.db.session.add(trans)

            
            indx += 1
        
        
    def Convert(self):
        
        self.ReadMindata()
        self.ReadTSdata()
        
    def getDatabase(self):
        
        return self.db
    

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('natoms', help='Number of atoms in system', type = int)
    parser.add_argument('--Database','-d', help = 'Name of database to write into', type = str)
    parser.add_argument('--Mindata','-m', help = 'Name of min.data file', type = str)
    parser.add_argument('--Tsdata','-t', help = 'Name of ts.data file', type = str)
    parser.add_argument('--Pointsmin','-p', help = 'Name of points.min file', type = str)
    parser.add_argument('--Pointsts','-q', help = 'Name of points.ts file', type = str)
    args = parser.parse_args()
    
    cv = Convert(args.Database)

    cv.initDatabase()
    cv.setAccuracy()
    
    cv.setNatoms(args.natoms)
     
    cv.setMindata(args.Mindata)
    cv.setTsdata(args.Tsdata)
    cv.setPointsts(args.Pointsmin)
    cv.setPointsmin(args.Pointsts)
     
    cv.Convert()
    cv.db.session.commit()


    

