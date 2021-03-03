"""
tools for reading and writing OPTIM input and output files
"""
from __future__ import print_function

import numpy as np
from pele.storage import Minimum, TransitionState

_id_count = 0

# class UnboundMinimum(object):
# """
# a class to duplicate some of the functionality of the Minimum class
# """
# energy = None
# coords = None
# fvib = None
# pgorder = None
# #    _id_count = 0
# def __init__(self, energy, coords):
# self.energy = energy
# self.coords = coords
# global _id_count
# self.id() = _id_count
# _id_count += 1
#
#    def __eq__(self, m):
#        """m can be integer or Minima object"""
#        assert self.id() is not None
#        if isinstance(m, UnboundMinimum):
#            assert m.id() is not None
#            return self.id() == m.id()
#        else:
#            return self.id() == m
#        
#    def __hash__(self):
#        assert self.id() is not None
#        return self.id()

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
        coords = np.fromfile(fin, dtype=np.dtype(endianness + "d"))
    if ndof is not None:
        if len(coords) % ndof != 0:
            raise Exception("number of double precision variables read from %s (%s) is not divisible by ndof (%d)" %
                            (fname, len(coords), ndof))
        #    print coords
    return coords.reshape(-1)

def write_points_min_ts(fout, x, endianness="="):
    """
    write coords to a points.min or a points.ts file
    """
    x = np.asarray(x.ravel(), dtype=np.dtype(endianness + "d"))
    x.tofile(fout)

class OptimDBConverter(object):
    """
    converts PATHSAMPLE to pele database

    Parameters
    ----------
    database : pele Database
        the minima and transition states will be place in here
    ndof : int, optional
        for testing to make sure the number of floats read is a multiple of ndof
    mindata, tsdata, pointsmin, pointsts : str
        the files to read from.  The files contain

            points.min : the coordinates of the minima in binary format
            min.data   : additional information about the minima (like the energy)
            points.ts  : the coordinates of the transition states
            min.ts     : additional information about transition states (like which minima they connect)

    endianness : str
        define the endianness of the binary data. can be "=", "<", ">"
    assert_coords : bool
        If this is True the conversion will abort if the coordinate conversion doesn't work.
        Set this to false if you only care about the minima and ts metadata.

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


    """

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

    def setAccuracy(self, accuracy=0.000001):
        self.db.accuracy = accuracy

    def ReadMinDataFast(self):
        """read min.data file
        
        this method uses bulk database inserts.  It is *MUCH* faster this way, but 
        you have to be careful that this and the Minimum object stays in sync.  e.g.
        minimum.invalid must be set to false manually here.
        """
        print("reading from", self.mindata)
        indx = 0
        #        f_len = file_len(self.mindata)
        minima_dicts = []
        for line in open(self.mindata, 'r'):
            sline = line.split()

            # get the coordinates corresponding to this minimum
            if self.pointsmin_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsmin_data[indx, :]


            # read data from the min.data line            
            e, fvib = list(map(float, sline[:2]))  # energy and vibrational free energy
            pg = int(sline[2])  # point group order

            # create the minimum object and attach the data
            # must add minima like this.  If you use db.addMinimum()
            # some minima with similar energy might be assumed to be duplicates
            min_dict = dict(energy=e, coords=coords, invalid=False,
                            fvib=fvib, pgorder=pg
            )
            minima_dicts.append(min_dict)

            indx += 1
        #            if indx % 50 == 0:
        #                self.db.session.commit()

        self.db.engine.execute(Minimum.__table__.insert(), minima_dicts)
        self.db.session.commit()

        print("--->finished loading %s minima" % indx)


    def ReadMindata(self):  # pragma: no cover
        print("reading from", self.mindata)
        indx = 0
        #        f_len = file_len(self.mindata)
        self.index2min = dict()
        for line in open(self.mindata, 'r'):
            sline = line.split()

            # get the coordinates corresponding to this minimum
            if self.pointsmin_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsmin_data[indx, :]

            # read data from the min.data line            
            e, fvib = list(map(float, sline[:2]))  # energy and vibrational free energy
            pg = int(sline[2])  # point group order

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

        print("--->finished loading %s minima" % indx)

    def ReadTSdataFast(self):
        """read ts.data file
        
        this method uses bulk database inserts.  It is *MUCH* faster this way, but 
        you have to be careful that this and the TransitionState object stays in sync.  e.g.
        ts.invalid must be set to false manually here.

        """
        print("reading from", self.tsdata)

        indx = 0
        ts_dicts = []
        for line in open(self.tsdata, 'r'):
            sline = line.split()

            # get the coordinates corresponding to this minimum
            if self.pointsts_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsts_data[indx, :]

            # read data from the min.ts line            
            e, fvib = list(map(float, sline[:2]))  # get energy and fvib
            pg = int(sline[2])  # point group order
            m1indx, m2indx = list(map(int, sline[3:5]))
            #            m1indx -= 1
            #            m2indx -= 1
            #            min1 = self.index2min[m1indx - 1] # minus 1 for fortran indexing
            #            min2 = self.index2min[m2indx - 1] # minus 1 for fortran indexing

            # must add transition states like this.  If you use db.addtransitionState()
            # some transition states might be assumed to be duplicates
            tsdict = dict(energy=e, coords=coords, invalid=False,
                          fvib=fvib, pgorder=pg,
                          _minimum1_id=m1indx,
                          _minimum2_id=m2indx
            )
            ts_dicts.append(tsdict)

            indx += 1
        #            if indx % 50 == 0:
        #                self.db.session.commit()
        self.db.engine.execute(TransitionState.__table__.insert(), ts_dicts)
        self.db.session.commit()

        print("--->finished loading %s transition states" % indx)


    def ReadTSdata(self):  # pragma: no cover
        print("reading from", self.tsdata)

        indx = 0
        for line in open(self.tsdata, 'r'):
            sline = line.split()

            # get the coordinates corresponding to this minimum
            if self.pointsts_data is None:
                coords = np.zeros(1)
            else:
                coords = self.pointsts_data[indx, :]

            # read data from the min.ts line            
            e, fvib = list(map(float, sline[:2]))  # get energy and fvib
            pg = int(sline[2])  # point group order
            m1indx, m2indx = list(map(int, sline[3:5]))

            min1 = self.index2min[m1indx - 1]  # minus 1 for fortran indexing
            min2 = self.index2min[m2indx - 1]  # minus 1 for fortran indexing

            # must add transition states like this.  If you use db.addtransitionState()
            # some transition states might be assumed to be duplicates
            trans = TransitionState(e, coords, min1, min2)

            trans.fvib = fvib
            trans.pgorder = pg

            indx += 1
            self.db.session.add(trans)
            if indx % 50 == 0:
                self.db.session.commit()

        print("--->finished loading %s transition states" % indx)

    def read_points_min(self):
        print("reading from", self.pointsmin)
        coords = read_points_min_ts(self.pointsmin, self.ndof, endianness=self.endianness)
        if coords.size == 0:
            raise Exception(self.pointsmin + " is empty")
        if self.ndof is None:
            # try to get the number of minima from the min.data file
            nminima = sum((1 for _ in open(self.mindata, "r")))
            assert len(coords.shape) == 1
            if coords.size % nminima != 0:
                raise ValueError("the number of data points in %s is not divisible by %s the number of minima in %s"
                                 % (self.mindata, coords.size, nminima))
            self.ndof = coords.size // nminima
            print("read %s minimum coordinates of length %s" % (nminima, self.ndof))
        self.pointsmin_data = coords.reshape([-1, self.ndof])

    def read_points_ts(self):
        print("reading from", self.pointsts)
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
        self.ReadMinDataFast()
        self.db.session.commit()

    def load_transition_states(self):
        try:
            self.read_points_ts()
        except IOError:
            if self.no_coords_ok:
                self.pointsts_data = None
            else:
                raise

        self.ReadTSdataFast()
        self.db.session.commit()

    def convert(self):
        """convert pathsample database to pele database"""
        self.load_minima()
        self.load_transition_states()

    def convert_no_coords(self):
        """convert pathsample database to pele database without loading coordinates"""
        self.pointsmin_data = None
        self.pointsts_data = None
        self.ReadMinDataFast()
        self.ReadTSdataFast()

class WritePathsampleDB(object):
    """
    converts PATHSAMPLE to pele database

    Parameters
    ----------
    database : pele Database
        the minima and transition states will be read from here
    mindata, tsdata, pointsmin, pointsts : str
        the files to writen to.  The files contain

            points.min : the coordinates of the minima in binary format
            min.data   : additional information about the minima (like the energy)
            points.ts  : the coordinates of the transition states
            min.ts     : additional information about transition states (like which minima they connect)

    endianness : str
        define the endianness of the binary data. can be "=", "<", ">"

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


    """

    def __init__(self, database, mindata="min.data",
                 tsdata="ts.data", pointsmin="points.min", pointsts="points.ts",
                 endianness="=", assert_coords=True):
        self.db = database
        self.mindata = mindata
        self.tsdata = tsdata
        self.pointsmin = pointsmin
        self.pointsts = pointsts
        self.endianness = endianness
    
    def write_min_data_ts_data(self):

        # write minima ordered by energy
        minima_labels = dict()
        import sqlalchemy.orm
        with open(self.pointsmin, "wb") as point_out: 
            with open(self.mindata, "w") as data_out:
                
                minima_iter = self.db.session.query(Minimum).\
                            options(sqlalchemy.orm.undefer("coords")).order_by(Minimum.energy)
                for label, m in enumerate(minima_iter):
                    minima_labels[m.id()] = label + 1 # +1 so it starts with 1
                    fvib = m.fvib
                    if fvib is None:
                        fvib = 1.
                    pgorder = m.pgorder
                    if pgorder is None:
                        pgorder = 1
                    
                    data_out.write("{} {} {} 1 1 1\n".format(m.energy,
                                                              fvib,
                                                              pgorder))
                    write_points_min_ts(point_out, m.coords, endianness=self.endianness)
        
        del m
        
        # write trasnition_states ordered by energy
        with open(self.pointsts, "wb") as point_out: 
            with open(self.tsdata, "w") as data_out:
                
                ts_iter = self.db.session.query(TransitionState).\
                            options(sqlalchemy.orm.undefer("coords"))
                for ts in ts_iter:
                    m1_label = minima_labels[ts._minimum1_id]
                    m2_label = minima_labels[ts._minimum2_id]

                    fvib = ts.fvib
                    if fvib is None:
                        fvib = 1.
                    pgorder = ts.pgorder
                    if pgorder is None:
                        pgorder = 1

                    data_out.write("{energy} {fvib} {pgorder} {min1} {min2} 1 1 1\n".format(
                        energy=ts.energy, fvib=fvib, pgorder=pgorder, 
                        min1=m1_label, min2=m2_label))
                    write_points_min_ts(point_out, ts.coords, endianness=self.endianness)
        
        
    
    def write_db(self):
        self.write_min_data_ts_data()
        
    

