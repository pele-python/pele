"""
tools for reading and writing OPTIM input and output files
"""

import numpy as np
from pele.storage import Database, Minimum, TransitionState

_id_count = 0

class UnboundMinimum(object):
    """
    a class to duplicate some of the functionality of the Minimum class
    """
    energy = None
    coords = None
    fvib = None
    pgorder = None
#    _id_count = 0
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = coords
        global _id_count
        self._id = _id_count
        _id_count += 1

    def __eq__(self, m):
        """m can be integer or Minima object"""
        assert self._id is not None
        if isinstance(m, UnboundMinimum):
            assert m._id is not None
            return self._id == m._id
        else:
            return self._id == m
        
    def __hash__(self):
        assert self._id is not None
        return self._id

def read_min_data(fname="min.data"):
    """
    return a list of minima with data read in from a min.data file
    
    Parameters
    ----------
    fname : str
        the file name of the min.data file
    
    Returns
    -------
    list of UnboundMinimum objects.
    """
    minima = []
    
    with open(fname, "r") as fin:
        for line in fin:
            sline = line.split()
            energy = float(sline[0])
            coords = np.array([0.])
            m = UnboundMinimum(energy, coords)
            
            m.fvib = float(sline[1])
            m.pgorder = int(sline[2])
            
            minima.append(m)
    return minima


def read_points_min_ts(fname, ndof=None, endianness="="):
    """
    read coords from a points.min or a points.ts file
    
    Notes
    -----
    the files were written with fortran code that looks something like this::
    
        NOPT = 3 * NATOMS
        INQUIRE(IOLENGTH=NDUMMY) COORDS(1:NOPT)
        OPEN(13,FILE='points.min,ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=NDUMMY)
        DO J1=1,NMIN
            WRITE(13,REC=J1) COORDS(1:NOPT)
        ENDDO
        CLOSE(13) 
    
    This means the data is stored without any header information.  
    It is just a long list of double precision floating point numbers.
    
    Note that some fortran compilers use different endiness for the data.  If
    the coordinates comes out garbage this is probably the problem.  The solution
    is to pass a different data type
    
    dtype=np.dtype("<d")  # for little-endian double precision
    dtype=np.dtype(">d")  # for big-endian double precision
    
    Parameters
    ----------
    fname : str
        filenname to read from
    ndof : int, optional
        for testing to make sure the number of floats read is a multiple of ndof
    endianness : str
        define the endianness of the data. can be "=", "<", ">" 
     
    """
    with open(fname, "rb") as fin:
        coords = np.fromfile(fin, dtype=np.dtype(endianness+"d"))
    if ndof is not None:
        if len(coords) % ndof != 0:
                raise Exception("number of double precision variables read from %s (%s) is not divisible by ndof (%d)" % 
                                (fname, len(coords), ndof) )
#    print coords
    return coords.reshape(-1)

class OptimDBConverter(object):
    '''
    Converts old OPTIM to pele database
    '''
    def __init__(self, database, ndof=None, mindata="min.data", 
                  tsdata="ts.data", pointsmin="points.min", pointsts="points.ts",
                  endianness="=", assert_coords=True):
        self.db = database
        self.ndof = ndof
        self.mindata = mindata
        self.tsdata = tsdata
        self.pointsmin = pointsmin
        self.pointsts = pointsts
        self.endianness = endianness
        self.no_coords_ok = not assert_coords

    def setAccuracy(self,accuracy = 0.000001):
        self.db.accuracy = accuracy
        
    def ReadMindata(self):
        print "reading from", self.mindata
        indx = 0
#        f_len = file_len(self.mindata)
        self.index2min = dict()
        for line in open(self.mindata,'r'):
            sline = line.split()
            
            # get the coordinates corresponding to this minimum
            if self.pointsmin_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsmin_data[indx,:]
            
            # read data from the min.data line            
            e, fvib = map(float,sline[:2]) # energy and vibrational free energy
            pg = int(sline[2]) # point group order
            
            # create the minimum object and attach the data
            # must add minima like this.  If you use db.addMinimum()
            # some minima with similar energy might be assumed to be duplicates
            min1 = Minimum(e, coords)

            min1.fvib = fvib
            min1.pgorder = pg
            
            self.index2min[indx] = min1

            indx += 1
            self.db.session.add(min1)
            if indx % 50 == 0:
                self.db.session.commit()

    def ReadTSdata(self):
        print "reading from", self.tsdata

        indx = 0
        for line in open(self.tsdata,'r'):
            sline = line.split()

            # get the coordinates corresponding to this minimum
            if self.pointsts_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsts_data[indx,:]

            # read data from the min.ts line            
            e, fvib = map(float, sline[:2]) # get energy and fvib
            pg = int(sline[2]) # point group order
            m1indx, m2indx = map(int, sline[3:5]) 
            
            min1 = self.index2min[m1indx - 1] # minus 1 for fortran indexing
            min2 = self.index2min[m2indx - 1] # minus 1 for fortran indexing

            # must add transition states like this.  If you use db.addtransitionState()
            # some transition states might be assumed to be duplicates
            trans = TransitionState(e, coords, min1, min2)
            
            trans.fvib = fvib
            trans.pgorder = pg
            
            indx += 1
            self.db.session.add(trans)
            if indx % 50 == 0:
                self.db.session.commit()
        
    def read_points_min(self):
        print "reading from", self.pointsmin
        coords = read_points_min_ts(self.pointsmin, self.ndof, endianness=self.endianness)
        if self.ndof is None:
            nminima = sum((1 for line in open(self.mindata, "r")))
            assert len(coords.shape) == 1
            if coords.size % nminima != 0:
                raise ValueError("the number of data points in %s is not divisible by %s the number of minima in %s"
                                 % (self.mindata, coords.size, nminima))
            self.ndof = coords.size / nminima
            print nminima, self.ndof, nminima*self.ndof, coords.size
        self.pointsmin_data = coords.reshape([-1, self.ndof])

    def read_points_ts(self):
        print "reading from", self.pointsts
        coords = read_points_min_ts(self.pointsts, self.ndof, endianness=self.endianness)
        self.pointsts_data = coords.reshape([-1, self.ndof])
                  
    def load_minima(self):
        try:
            self.read_points_min()
        except IOError:
            if self.no_coords_ok:
                self.pointsmin_data = None
            else:
                raise
        self.ReadMindata()
        self.db.session.commit()
    
    def load_transition_states(self):
        try:
            self.read_points_ts()
        except IOError:
            if self.no_coords_ok:
                self.pointsts_data = None
            else:
                raise

        self.ReadTSdata()
        self.db.session.commit()
    
    def convert(self):
        
        self.load_minima()
        self.load_transition_states()
