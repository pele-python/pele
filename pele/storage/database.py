"""Database for simulation data in a relational database
"""
from __future__ import print_function
import threading
import os

import numpy as np

from sqlalchemy import create_engine, and_, or_
from sqlalchemy.orm import sessionmaker, undefer
from sqlalchemy import Column, Integer, Float, PickleType, String
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, deferred
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.schema import Index

from pele.utils.events import Signal

__all__ = ["Minimum", "TransitionState", "Database"]

_schema_version = 2
verbose=False

Base = declarative_base()

# the following taken from http://www.rfk.id.au/blog/entry/preparing-pyenchant-for-python-3/
try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = basestring


class Minimum(Base):
    """
    The Minimum class represents a minimum in the database.

    Parameters
    ----------
    energy : float
    coords : numpy array
        coordinates

    Attributes
    ----------
    energy :
        the energy of the minimum
    coords :
        the coordinates of the minimum.  This is stored as a pickled numpy
        array which SQL interprets as a BLOB.
    fvib :
        the log product of the squared normal mode frequencies.  This is used in
        the free energy calcualations
    pgorder :
        point group order
    invalid :
        a flag that can be used to indicate a problem with the minimum.  E.g. if
        the Hessian has more zero eigenvalues than expected.
    user_data :
        Space to store anything that the user wants.  This is stored in SQL
        as a BLOB, so you can put anything here you want as long as it's serializable.
        Usually a dictionary works best.


    Notes
    -----
    To avoid any double entries of minima and be able to compare them,
    only use `Database.addMinimum()` to create a minimum object.

    See Also
    --------
    Database, TransitionState
    """
    __tablename__ = 'tbl_minima'

    _id = Column(Integer, primary_key=True)
    energy = Column(Float) 
    # deferred means the object is loaded on demand, that saves some time / memory for huge graphs
    coords = deferred(Column(PickleType))
    '''coordinates of the minimum'''
    fvib = Column(Float)
    """log product of the squared normal mode frequencies"""
    pgorder = Column(Integer)
    """point group order"""
    invalid = Column(Integer)
    """flag indicating if the minimum is invalid"""
    user_data = deferred(Column(PickleType))
    """this can be used to store information about the minimum"""
    
    
    def __init__(self, energy, coords):
        self.energy = energy
        self.coords = np.copy(coords)
        self.invalid = False

    def id(self):
        """return the sql id of the object"""
        return self._id

    def __eq__(self, m):
        """m can be integer or Minima object"""
        assert self.id() is not None
        if isinstance(m, Minimum):
            assert m.id() is not None
            return self.id() == m.id()
        else:
            return self.id() == m
        
    def __hash__(self):
        _id = self.id()
        assert _id is not None
        return _id

#    transition_states = relationship("transition_states", order_by="transition_states.id", backref="minima")
    
class TransitionState(Base):
    """Transition state object

    The TransitionState class represents a saddle point in the database.

    Parameters
    ----------
    energy : float
    coords : numpy array
    min1 : Minimum object
        first minimum
    min2 : Minimum object
        first minimum
    eigenval : float, optional
        lowest (single negative) eigenvalue of the saddle point
    eigenvec : numpy array, optional
        eigenvector which corresponds to the negative eigenvalue
    fvib : float
        log product of squared frequencies for free energy calculation
    pgorder : integer
        point group order



    Attributes
    ----------
    energy :
        The energy of the transition state
    coords :
        The coordinates of the transition state.  This is stored as a pickled numpy
        array which SQL interprets as a BLOB.
    fvib :
        The log product of the squared normal mode frequencies.  This is used in
        the free energy calcualations
    pgorder :
        The point group order
    invalid :
        A flag that is used to indicate a problem with the transition state.  E.g. if
        the Hessian has more than one negaive eigenvalue then it is a higher order saddle.
    user_data :
        Space to store anything that the user wants.  This is stored in SQL
        as a BLOB, so you can put anything here you want as long as it's serializable.
        Usually a dictionary works best.
    minimum1, minimum2 :
        These returns the minima on either side of the transition state
    eigenvec :
        The vector which points along the direction crossing the transition state.
        This is the eigenvector of the lowest non-zero eigenvalue.
    eigenval :
        The eigenvalue corresponding to `eigenvec`.  A.k.a. the curvature
        along the direction given by `eigenvec`

    Notes
    -----
    To avoid any double entries and be able to compare them, only use
    Database.addTransitionState to create a TransitionStateobject.

    programming note: The functions in the database require that
    ts.minimum1._id < ts.minimum2._id.  This will be handled automatically
    by the database, but we must remember not to screw it up

    See Also
    --------
    Database, Minimum
    """
    __tablename__ = "tbl_transition_states"
    _id = Column(Integer, primary_key=True)
    
    energy = Column(Float)
    '''energy of transition state'''
    
    coords = deferred(Column(PickleType))
    '''coordinates of transition state'''
    
    _minimum1_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum1 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum1_id")
    '''first minimum which connects to transition state'''
    
    _minimum2_id = Column(Integer, ForeignKey('tbl_minima._id'))
    minimum2 = relationship("Minimum",
                            primaryjoin="Minimum._id==TransitionState._minimum2_id")
    '''second minimum which connects to transition state'''
    
    eigenval = Column(Float)
    '''coordinates of transition state'''

    eigenvec = deferred(Column(PickleType))
    '''coordinates of transition state'''

    fvib = Column(Float)
    """log product of the squared normal mode frequencies"""
    pgorder = Column(Integer)
    """point group order"""
    invalid = Column(Integer)
    """flag indicating if the transition state is invalid"""
    user_data = deferred(Column(PickleType))
    """this can be used to store information about the transition state """


    def __init__(self, energy, coords, min1, min2, eigenval=None, eigenvec=None):
        assert min1.id() is not None
        assert min2.id() is not None
        
        self.energy = energy
        self.coords = np.copy(coords)
        if min1.id() < min2.id():
            self.minimum1 = min1
            self.minimum2 = min2
        else:
            self.minimum1 = min2
            self.minimum2 = min1
        
        if eigenvec is not None:
            self.eigenvec = np.copy(eigenvec)
        self.eigenval = eigenval
        self.invalid = False

    def id(self):
        """return the sql id of the object"""
        return self._id



class SystemProperty(Base):
    """table to hold system properties like potential parameters and number of atoms
    
    The properties can be stored as integers, floats, strings, or a pickled object.
    Only one of the property value types should be set for each property.
    """
    __tablename__ = "tbl_system_property"
    _id = Column(Integer, primary_key=True)

    property_name = Column(String)
    int_value = Column(Integer)
    float_value = Column(Float)
    string_value = Column(String)
    pickle_value = deferred(Column(PickleType))

    def __init__(self, property_name):
        self.property_name = property_name
        
    def name(self):
        return self.property_name
    
    def _values(self):
        """return a dictionary of the values that are not None"""
        values = dict(int_value=self.int_value, float_value=self.float_value, 
                      string_value=self.string_value, pickle_value=self.pickle_value)
        values = dict([(k,v) for k,v in values.items() if v is not None])
        return values
    
    def value(self):
        """return the property value"""
        actual_values = [v for v in self._values().values() if v is not None]
        if len(actual_values) == 1:
            return actual_values[0]
        elif len(actual_values) == 0:
            return None
        elif len(actual_values) > 1:
            print("SystemProperty: multiple property values are set")
            return actual_values
        return None
    
    def item(self):
        """return a tuple of (name, value)"""
        return self.name(), self.value()
            

Index('idx_transition_states', TransitionState.__table__.c._minimum1_id, TransitionState.__table__.c._minimum2_id)
Index('idx_minimum_energy', Minimum.__table__.c.energy)
Index('idx_transition_state_energy', Minimum.__table__.c.energy)


class MinimumAdder(object):
    """This class manages adding minima to the database
    
    Parameters
    ----------
    db : database object
    Ecut: float, optional
        energy cutoff, don't add minima which are higher in energy
    max_n_minima : int, optional
        keep only the max_n_minima with the lowest energies. If E is greater
        than the minimum with the highest energy in the database, then don't add
        this minimum and return None.  Else add this minimum and delete the minimum
        with the highest energy.  if max_n_minima < 0 then it is ignored.
    commit_interval : int, optional
        Commit the database changes to the hard drive every `commit_interval` steps.
        Committing too frequently can be slow, this is used to speed things up.
    """
    def __init__(self, db, Ecut=None, max_n_minima=None, commit_interval=1):
        self.db = db
        self.Ecut = Ecut
        self.max_n_minima = max_n_minima
        self.commit_interval = commit_interval
        self.count = 0

    def __call__(self, E, coords):
        """this is called to add a minimum to the database"""
        if self.Ecut is not None:
            if E > self.Ecut:
                return None
        commit = self.count % self.commit_interval == 0
        self.count += 1
        return self.db.addMinimum(E, coords, max_n_minima=self.max_n_minima, 
                                  commit=commit)
    
    def __del__(self):
        """ensure that all the changes to the database are committed to the hard drive
        """
        if self.commit_interval != 1:
            self.db.session.commit()

def _compare_properties(prop, v2):
    v1 = prop.value()
    try:
        return bool(v1 == v2)
    except Exception:
        pass
    
    try:
        # see if they are numpy arrays
        return np.all(v1 == v2)
    except:
        pass
    
    print("warning, could not compare value", v2, "with", v1)
    return False
        
    
        
            
        
        
        

class Database(object):
    """Database storage class

    The Database class handles the connection to the database. It has functions to create new Minima and
    TransitionState objects. The objects are persistent in the database and exist as
    soon as the Database class in connected to the database. If any value in the objects is changed,
    the changes are automatically persistent in the database (TODO: be careful, check commit transactions, ...)

    Database uses SQLAlchemy to connect to the database. Check the web page for available connectors. Unless
    you know better, the standard sqlite should be used. The database can be generated in memory (default) or
    written to a file if db is specified when creating the class.

    Parameters
    ----------
    db : string, optional
        filename of new or existing database to connect to.  default creates
        new database in memory.
    accuracy : float, optional
        energy tolerance to count minima as equal
    connect_string : string, optional
        connection string, default is sqlite database
    compareMinima : callable, `bool = compareMinima(min1, min2)`, optional
        called to determine if two minima are identical.  Only called
        if the energies are within `accuracy` of each other.
    createdb : boolean, optional
        create database if not exists, default is true

    Attributes
    ----------
    engine : sqlalchemy database engine
    session : sqlalchemy session

    accuracy : float
    on_minimum_removed : signal
        called when a minimum is removed from the database
    on_minimum_added : signal
        called when a new, unique, minimum is added to the database
    on_ts_removed : signal
        called when a transition_state is removed from the database
    on_ts_added : signal
        called when a new, unique, transition state is added to the database
    compareMinima

    Examples
    --------

    >>> from pele.storage import Database
    >>> db = Database(db="test.db")
    >>> for energy in np.random.random(10):
    >>>     a.addMinimum(energy, np.random.random(10))
    >>>
    >>> for minimum in database.minima():
    >>>     print minimum.energy

    See Also
    --------
    Minimum
    TransitionState

    """
    engine = None
    session = None
    connection = None
    accuracy = 1e-3
    compareMinima=None
        
    def __init__(self, db=":memory:", accuracy=1e-3, connect_string='sqlite:///%s',
                 compareMinima=None, createdb=True):
        self.accuracy=accuracy
        self.compareMinima = compareMinima

        if not os.path.isfile(db) or db == ":memory:":
            newfile = True
            if not createdb:
                raise IOError("createdb is False, but database does not exist (%s)" % db)
        else:
            newfile = False

        # set up the engine which will manage the backend connection to the database
        self.engine = create_engine(connect_string % db, echo=verbose)

        if not newfile and not self._is_pele_database():
            raise IOError("existing file (%s) is not a pele database." % db)
        

        # set up the tables and check the schema version
        if newfile:
            self._set_schema_version()
        self._check_schema_version()
        self._update_schema()
#         self._check_schema_version_and_create_tables(newfile)

        # set up the session which will manage the frontend connection to the database
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        
        # these functions will be called when a minimum or transition state is 
        # added or removed
        self.on_minimum_added = Signal()
        self.on_minimum_removed = Signal()
        self.on_ts_added = Signal()
        self.on_ts_removed = Signal()
        
        self.lock = threading.Lock()
        self.connection = self.engine.connect()

    def _is_pele_database(self):
        conn = self.engine.connect()
        result = True
        if (not self.engine.has_table("tbl_minima") or
                not self.engine.has_table("tbl_transition_states")):
            result = False
        conn.close()
        return result

    def _set_schema_version(self):
        conn = self.engine.connect()
        conn.execute("PRAGMA user_version = %d;"%_schema_version)
        conn.close()

    def _update_schema(self):
        conn = self.engine.connect()
        Base.metadata.create_all(bind=self.engine)
        conn.close()

    def _check_schema_version(self):
        conn = self.engine.connect()
        result=conn.execute("PRAGMA user_version;")
        schema = result.fetchone()[0]
        result.close()
        conn.close()
        if _schema_version != schema:
            raise IOError("database schema outdated, current (newest) version: "
                          "%d (%d). Please use migrate_db.py in pele/scripts to update database"%(schema, _schema_version))


    def _highest_energy_minimum(self):
        """return the minimum with the highest energy"""
        candidates = self.session.query(Minimum).order_by(Minimum.energy.desc()).\
            limit(1).all()
        return candidates[0]
    
    def get_lowest_energy_minimum(self):
        """return the minimum with the lowest energy"""
        candidates = self.session.query(Minimum).order_by(Minimum.energy).\
            limit(1).all()
        return candidates[0]
    
    def findMinimum(self, E, coords):
        candidates = self.session.query(Minimum).\
            options(undefer("coords")).\
            filter(Minimum.energy > E-self.accuracy).\
            filter(Minimum.energy < E+self.accuracy)
        
        new = Minimum(E, coords)
        
        for m in candidates:
            if self.compareMinima:
                if not self.compareMinima(new, m):
                    continue
            return m
        return None
        
    def addMinimum(self, E, coords, commit=True, max_n_minima=-1, pgorder=None, fvib=None):
        """add a new minimum to database
        
        Parameters
        ----------
        E : float
        coords : numpy.array
            coordinates of the minimum
        commit : bool, optional
            commit changes to database
        max_n_minima : int, optional
            keep only the max_n_minima with the lowest energies. If E is greater
            than the minimum with the highest energy in the database, then don't add
            this minimum and return None.  Else add this minimum and delete the minimum
            with the highest energy.  if max_n_minima < 0 then it is ignored.

        Returns
        -------
        minimum : Minimum
            minimum which was added (not necessarily a new minimum)
            
        """
        self.lock.acquire()
        # undefer coords because it is likely to be used by compareMinima and
        # it is slow to load them individually by accessing the database repetitively.
        candidates = self.session.query(Minimum).\
            options(undefer("coords")).\
            filter(Minimum.energy.between(E-self.accuracy, E+self.accuracy))
        
        new = Minimum(E, coords)
        
        for m in candidates:
            if self.compareMinima:
                if not self.compareMinima(new, m):
                    continue
            self.lock.release()
            return m

        if max_n_minima is not None and max_n_minima > 0:
            if self.number_of_minima() >= max_n_minima:
                mmax = self._highest_energy_minimum()
                if E >= mmax.energy:
                    # don't add the minimum
                    self.lock.release() 
                    return None
                else:
                    # remove the minimum with the highest energy and continue
                    self.removeMinimum(mmax, commit=commit)

        if fvib is not None:
            new.fvib = fvib
        if pgorder is not None:
            new.pgorder = pgorder
        self.session.add(new)
        if commit:
            self.session.commit()
        
        self.lock.release()
        
        self.on_minimum_added(new)
        return new
        
    def getMinimum(self, mid):
        """return the minimum with a given id"""
        return self.session.query(Minimum).get(mid)
        
    def addTransitionState(self, energy, coords, min1, min2, commit=True, 
                           eigenval=None, eigenvec=None, pgorder=None, fvib=None):
        """Add transition state object
        
        Parameters
        ----------
        energy : float
            energy of transition state
        coords : numpy array
            coordinates of transition state
        min1, min2 : Minimum
            minima on either side of the transition states
        eigenval : float
            the eigenvalue (curvature) across the transition state
        eigenvec : numpy array
            the eigenvector associated with eigenval
        commit : bool
            commit changes to sql database
        
        Returns
        -------
        ts : TransitionState
            the transition state object (not necessarily new)
        """
        m1, m2 = min1, min2
        if m1.id() > m2.id():
            m1, m2 = m2, m1
        candidates = self.session.query(TransitionState).\
            options(undefer("coords")).\
            filter(or_(
                       and_(TransitionState.minimum1==m1, 
                            TransitionState.minimum2==m2),
                       and_(TransitionState.minimum1==m2, 
                            TransitionState.minimum2==m1),
                       )).\
            filter(TransitionState.energy.between(energy-self.accuracy,  energy+self.accuracy))
        
        for m in candidates:
            return m

        new = TransitionState(energy, coords, m1, m2, eigenval=eigenval, eigenvec=eigenvec)
        
        if fvib is not None:
            new.fvib = fvib
        if pgorder is not None:
            new.pgorder = pgorder 
        self.session.add(new)
        if commit:
            self.session.commit()
        self.on_ts_added(new)
        return new

    def getTransitionState(self, min1, min2):
        """return the TransitionState between two minima
        
        Returns
        -------
        ts : None or TransitionState
        """
        m1, m2 = min1, min2
        candidates = self.session.query(TransitionState).\
            filter(or_(
                       and_(TransitionState.minimum1==m1, 
                            TransitionState.minimum2==m2),
                       and_(TransitionState.minimum1==m2, 
                            TransitionState.minimum2==m1),
                       ))

        for m in candidates:
            return m
        return None

    def getTransitionStatesMinimum(self, min1):
        """return all transition states connected to a minimum
        
        Returns
        -------
        ts : None or TransitionState
        """
        candidates = self.session.query(TransitionState).\
            filter(or_(TransitionState.minimum1==min1, 
                       TransitionState.minimum2==min1))

        return candidates.all()

    def getTransitionStateFromID(self, id_):
        """return the transition state with id id_"""
        return self.session.query(TransitionState).get(id_)
    
    def minima(self, order_energy=True):
        """return an iterator over all minima in database

        Parameters
        ----------
        order_energy : bool
            order the minima by energy

        Notes
        -----
        Minimum.coords is deferred in database queries.  If you need to access
        coords for multiple minima it is *much* faster to `undefer` before
        executing the query by, e.g.
        `session.query(Minimum).options(undefer("coords"))`
        """
        if order_energy:
            return self.session.query(Minimum).order_by(Minimum.energy).all()
        else:
            return self.session.query(Minimum).all()
    
    def transition_states(self, order_energy=False):
        """return an iterator over all transition states in database
        """
        if order_energy:
            return self.session.query(TransitionState).order_by(TransitionState.energy).all()
        else:
            return self.session.query(TransitionState).all()
    
    def minimum_adder(self, Ecut=None, max_n_minima=None, commit_interval=1):
        """wrapper class to add minima

        Since pickle cannot handle pointer to member functions, this class wraps the call to
        add minimum.

        Parameters
        ----------
        Ecut: float, optional
             energy cutoff, don't add minima which are higher in energy
        max_n_minima : int, optional
            keep only the max_n_minima with the lowest energies. If E is greater
            than the minimum with the highest energy in the database, then don't add
            this minimum and return None.  Else add this minimum and delete the minimum
            with the highest energy.  if max_n_minima < 0 then it is ignored.

        Returns
        -------
        handler: minimum_adder class
            minimum handler to add minima


        """
        return MinimumAdder(self, Ecut=Ecut, max_n_minima=max_n_minima,
                            commit_interval=commit_interval)
    
    def removeMinimum(self, m, commit=True):
        """remove a minimum from the database
        
        Remove a minimum and any objects (TransitionState) 
        pointing to that minimum.
        """
        # delete any transition states objects pointing to min2
        candidates = self.session.query(TransitionState).\
            filter(or_(TransitionState.minimum1 == m, 
                       TransitionState.minimum2 == m))
        candidates = list(candidates)
        for ts in candidates:
            self.on_ts_removed(ts)
            self.session.delete(ts)
        
        self.on_minimum_removed(m)
        # delete the minimum
        self.session.delete(m)
        if commit:
            self.session.commit()

    
    def mergeMinima(self, min1, min2):
        """merge two minima in the database
        
        min2 will be deleted and everything that 
        points to min2 will point to min1 instead.
        """
        # find all transition states for which ts.minimum1 is min2
        candidates = self.session.query(TransitionState).\
            filter(TransitionState.minimum1 == min2) 
        for ts in candidates:
            # should we check if this will duplicate an existing transition state?
            ts.minimum1 = min1
            if ts.minimum1.id() > ts.minimum2.id():
                ts.minimum1, ts.minimum2 = ts.minimum2, ts.minimum1
        
        # find all transition states for which ts.minimum2 is min2
        candidates = self.session.query(TransitionState).\
            filter(TransitionState.minimum2 == min2) 
        for ts in candidates:
            # should we check if this will duplicate an existing transition state?
            ts.minimum2 = min1
            if ts.minimum1.id() > ts.minimum2.id():
                ts.minimum1, ts.minimum2 = ts.minimum2, ts.minimum1
        
        self.session.delete(min2)
        self.session.commit()

    def remove_transition_state(self, ts, commit=True):
        """remove a transition states from the database
        """
        self.on_ts_removed(ts)
        self.session.delete(ts)
        if commit:
            self.session.commit()

    def number_of_minima(self):
        """return the number of minima in the database
        
        Notes
        -----
        This is much faster than len(database.minima()), but is is not instantaneous.  
        It takes a longer time for larger databases.  The first call to number_of_minima() 
        can be much faster than subsequent calls.  
        """
        return self.session.query(Minimum).count()

    def number_of_transition_states(self):
        """return the number of transition states in the database
        
        Notes
        -----
        see notes for number_of_minima()
        
        See Also
        --------
        number_of_minima
        """
        return self.session.query(TransitionState).count()

    def get_property(self, property_name):
        """return the minimum with a given name"""
        candidates = self.session.query(SystemProperty).\
            filter(SystemProperty.property_name == property_name)
        return candidates.first()

    def properties(self, as_dict=False):
        query = self.session.query(SystemProperty)
        if as_dict:
            return dict([p.item() for p in query])
        else:
            return query.all()

    def add_property(self, name, value, dtype=None, commit=True, overwrite=True):
        """add a system property to the database
        
        Parameters
        ----------
        name : string
            the name of the property
        value : object
            the value of the property
        dtype : string
            the datatype of the property.  This can be "int", "float", 
            "string", "pickle", or None.  If None, the datatype will be
            automatically determined.
        
        This could anything, such as a potential parameter, the number of atoms, or the
        list of frozen atoms. The properties can be stored as integers, floats, 
        strings, or a pickled object.  Only one of the property value types 
        should be set for each property.
        
        For a value of type "pickle", pass the object you want pickled, not 
        the pickled object.  We will do the pickling and unpickling for you.

        """
        new = self.get_property(name)
        if new is None:
            new = SystemProperty(name)
        else:
            # the database already has a property with this name, Try to determine if they are the same             
            same = _compare_properties(new, value)
            if not same:
                if not overwrite:
                    raise RuntimeError("property %s already exists and the value %s does not compare equal to the new value." % (new.item(), value))
                print("warning: overwriting old property", new.item())

        if dtype is None:
            # try to determine type of the value
            if isinstance(value, int):
                dtype = "int"
            elif isinstance(value, float):
                dtype = "float"
            elif isinstance(value, basestring):
                dtype = "string"
            else:
                dtype = "pickle"
        
        if dtype == "string":
            new.string_value = value
        elif dtype == "int":
            new.int_value = value
        elif dtype == "float":
            new.float_value = value
        elif dtype == "pickle":
            new.pickle_value = value
        else:
            raise ValueError('dtype must be one of "int", "float", "string", "pickle", or None')
            
        self.session.add(new)
        if commit:
            self.session.commit()
        return new
    
    def add_properties(self, properties, overwrite=True):
        """add multiple properties from a dictionary
        
        properties : dict
            a dictionary of (name, value) pairs.  The data type of the value
            will be determined automatically
        """
        for name, value in properties.items():
            self.add_property(name, value, commit=True, overwrite=overwrite)
            


def test_fast_insert(): # pragma: no cover
    """bulk inserts are *really* slow, we should add something along the lines of this
    answer to speed things up where needed
    
    http://stackoverflow.com/questions/11769366/why-is-sqlalchemy-insert-with-sqlite-25-times-slower-than-using-sqlite3-directly
    """
    db = Database()
    print(Minimum.__table__.insert())
    db.engine.execute(
                      Minimum.__table__.insert(),
                      [dict(energy=.01, coords=np.array([0.,1.]), invalid=False),
                      dict(energy=.02, coords=np.array([0.,2.]), invalid=False),
                      ]
                      )
    m1, m2 = db.minima()[:2]
    db.engine.execute(TransitionState.__table__.insert(),
                      [dict(energy=1., coords=np.array([1,1.]), _minimum1_id=m1.id(), 
                            _minimum2_id=m2.id())
                       ]
                      )
    for m in db.minima():
        print(m.id())
        print(m.energy)
        print(m.coords)
        print(m.invalid, bool(m.invalid))
        
    ts = db.transition_states()[0]
    print(ts.minimum1.energy)
    print(ts.minimum2.energy)
    print(ts.id())

if __name__ == "__main__":
    test_fast_insert()

