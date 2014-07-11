# distutils: language = c++
# distutils: sources = modified_fire.cpp
import numpy as np

from pele.potentials import _pele, _pythonpotential

#cimport numpy as np
#cimport cython
#from libcpp cimport bool as cbool

cdef class _Cdef_MODIFIED_FIRE_CPP(_pele_opt.GradientOptimizer):
    """This class is the python interface for the c++ MODIFIED_FIRE implementation
    """  
    cdef _pele.BasePotential pot
    
    def __cinit__(self, x0, potential, double dtstart = 0.1, double dtmax = 1, double maxstep=0.5, size_t Nmin=5, double finc=1.1, 
                  double fdec=0.5, double fa=0.99, double astart=0.1, double tol=1e-3, cbool stepback = True, 
                  int iprint=-1, energy=None, gradient=None, int nsteps=10000, int verbosity=0, events = None):
        
        if not issubclass(potential.__class__, _pele.BasePotential):
            if verbosity > 0:
                print "MODIFIED_FIRE_CPP: potential is not subclass of BasePotential; wrapping it.", potential
#               print "Wrapping the potential like this is dangerous.  All python exceptions will be ignored"
            potential = _pythonpotential.CppPotentialWrapper(potential)
        
        cdef _pele.BasePotential pot = potential
        cdef np.ndarray[double, ndim=1] x0c = np.array(x0, dtype=float)
        self.thisptr = shared_ptr[_pele_opt.cGradientOptimizer]( <_pele_opt.cGradientOptimizer*>
                        new cppMODIFIED_FIRE(pot.thisptr, _pele.Array[double](<double*> x0c.data, x0c.size),
                                             dtstart, dtmax, maxstep, Nmin, finc, fdec, fa, astart, tol, stepback) )
        
        self.thisptr.get().set_max_iter(nsteps)
        self.thisptr.get().set_verbosity(verbosity)
        self.thisptr.get().set_iprint(iprint)
        self.pot = pot
        
        cdef np.ndarray[double, ndim=1] g_  
        if energy is not None and gradient is not None:
            g_ = gradient
            self.thisptr.get().set_func_gradient(energy, _pele.Array[double](<double*> g_.data, g_.size))

        self.events = events
        if self.events is None: 
            self.events = []


class ModifiedFireCPP(_Cdef_MODIFIED_FIRE_CPP):
    """This class is the python interface for the c++ MODIFED_FIRE implementation.
    """
    
#     def reset(self, coords):
#         """do one iteration"""
#         _Cdef_MODIFIED_FIRE_CPP.reset(self, coords)
       
