import numpy as np
cimport numpy as np
from pele.potentials import _pele, _pythonpotential
from pele.potentials cimport _pele
from pele.optimize import Result
cimport cython
from pele.potentials import _pythonpotential
import sys
from libcpp cimport bool as cbool

# import the externally defined modified_fire implementation
cdef extern from "pele/modified_fire.h" namespace "pele":
    cdef cppclass cppMODIFIED_FIRE "pele::MODIFIED_FIRE":
        cppMODIFIED_FIRE(_pele.cBasePotential * , _pele.Array[double], 
                         double, double, double, size_t , double, double, 
                         double, double, double, cbool) except +
        void run() except +
        void run(int niter) except +
        void one_iteration() except +
        int stop_criterion_satisfied() except +
        void set_func_gradient(double energy, _pele.Array[double] grad) except +

        double get_f() except +
        _pele.Array[double] get_x() except +
        _pele.Array[double] get_g() except +
        
        void set_tol(double) except +
        void set_maxstep(double) except +
        void set_max_iter(int) except +
        void set_iprint(int) except +
        void set_verbosity(int) except +

        double get_rms() except +
        int get_nfev() except +
        int get_niter() except +
        int get_maxiter() except +
        int success() except +


cdef class _Cdef_MODIFIED_FIRE_CPP(object):
    """This class is the python interface for the c++ MODIFIED_FIRE implementation
    """
    cdef cppMODIFIED_FIRE *thisptr
    cdef _pele.BasePotential pot # this is stored so that the memory is not freed
    
    def __cinit__(self, x0, potential, double dtstart = 0.1, double dtmax = 1, double maxstep=0.5, size_t Nmin=5, double finc=1.1, 
                  double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-3, cbool stepback = True, 
                  int iprint=-1, energy=None, gradient=None, int nsteps=100000, int verbosity=0, events = None):
        
        if not issubclass(potential.__class__, _pele.BasePotential):
            if verbosity > 0:
                print "MODIFIED_FIRE_CPP: potential is not subclass of BasePotential; wrapping it.", potential
#               print "Wrapping the potential like this is dangerous.  All python exceptions will be ignored"
            potential = _pythonpotential.CppPotentialWrapper(potential)
        
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = <cppMODIFIED_FIRE*>new cppMODIFIED_FIRE(pot.thisptr, 
                                                               _pele.Array[double](<double*> x0c.data, x0c.size),
                                                               dtstart, dtmax, maxstep, Nmin, finc, fdec, fa, astart, tol, stepback)
        opt = self.thisptr
        opt.set_max_iter(nsteps)
        opt.set_verbosity(verbosity)
        opt.set_iprint(iprint)
        self.pot = pot # so that the memory is not freed
        
        cdef np.ndarray[double, ndim=1] g_  
        if energy is not None and gradient is not None:
            g_ = gradient
            self.thisptr.set_func_gradient(energy, _pele.Array[double](<double*> g_.data, g_.size))
        
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
        
    def run(self, niter=None):
        if niter is None:
            self.thisptr.run()
        else:
            self.thisptr.run(niter)
        return self.get_result()

    @cython.boundscheck(False)
    def get_result(self):
        """return a results object"""
        res = Result()
        
        cdef _pele.Array[double] xi = self.thisptr.get_x()
        cdef _pele.Array[double] gi = self.thisptr.get_g()
        cdef double *xdata = xi.data()
        cdef double *gdata = gi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        cdef np.ndarray[double, ndim=1, mode="c"] g = np.zeros(xi.size())
        cdef size_t i
        for i in xrange(xi.size()):
            x[i] = xdata[i]
            g[i] = gdata[i]

        #jake> it's anoying having to copy the data manually like this.
        # We can possibly use np.frombuffer(), thought I haven't gotten it to work.
        # We can also maybe use the c function PyArray_SimpleNewFromData.
        # In the meantime the loop won't be too slow if we use cython properly to speed it up.
        #cdef np.ndarray[double, ndim=1, mode="c"] g2
        #g2 = np.frombuffer(xi.data(), dtype=np.float64, count=xi.size())
        
        res.energy = self.thisptr.get_f()
        res.coords = x
        res.grad = g
        
        res.rms = self.thisptr.get_rms()        
        res.nsteps = self.thisptr.get_niter()
        res.nfev = self.thisptr.get_nfev()
        res.success = bool(self.thisptr.success())
        return res
    
    def one_iteration(self):
        self.thisptr.one_iteration()
    
    def stop_criterion_satisfied(self):
        return bool(self.thisptr.stop_criterion_satisfied())

    def get_maxiter(self):
        return self.thisptr.get_maxiter()

    def get_niter(self):
        return self.thisptr.get_niter()

class ModifiedFireCPP(_Cdef_MODIFIED_FIRE_CPP):
    """This class is the python interface for the c++ MODIFED_FIRE implementation.
    """
    def __init__(self, x0, potential, double dtstart=0.1, double dtmax=1, double maxstep=0.5, size_t Nmin=5,
                 double finc=1.1, double fdec=0.5, double fa=0.99, double astart=0.1, 
                 double tol=1e-3, cbool stepback = True, int iprint=-1, energy=None, gradient=None, 
                 int nsteps=100000, int verbosity=0,events=None):
        self._need_python = events is not None

        self.events = events
        if self.events is None: 
            self.events = []
        
    def one_iteration(self):
        """do one iteration"""
        _Cdef_MODIFIED_FIRE_CPP.one_iteration(self)
        res = self.get_result()
        for event in self.events:
            event(coords=res.coords, energy=res.energy, rms=res.rms)
    
    def run(self, niter=None):
        """run the lbfgs algorithm
        
        this overloads the underlying implementation so that the purely python 
        events can be called. If events are not called this simply calls the 
        purely cpp version.
        """
        if not self._need_python:
            return _Cdef_MODIFIED_FIRE_CPP.run(self, niter)

        if niter is None:
            niter = self.get_maxiter() - self.get_niter()

        for i in xrange(niter):
            if self.stop_criterion_satisfied():
                break
            self.one_iteration()
        
        return self.get_result()