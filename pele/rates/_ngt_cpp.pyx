"""
# distutils: language = C++
"""
import time

from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.list cimport list as stdlist

cdef extern from "pele/graph.hpp" namespace "pele":
    ctypedef unsigned long node_id

ctypedef pair[node_id, node_id] pair_t
ctypedef map[pair_t, double] rate_map_t

cdef extern from "pele/ngt.hpp" namespace "pele":
    cdef cppclass cNGT "pele::NGT":
        cNGT(rate_map_t &, stdlist[node_id] &, stdlist[node_id] &) except +
        void compute_rates() except +
        void compute_rates_and_committors() except +
        double get_rate_AB() except +
        double get_rate_BA() except +
        double get_rate_AB_SS() except +
        double get_rate_BA_SS() except +
        void set_node_occupation_probabilities(map[node_id, double] &) except +
        void set_debug() except +
        map[node_id, double] get_committors() except + # as reference ?




cdef class BaseNGT(object):
    """
    class to apply the graph reduction method for finding transition rates between two groups of nodes
    
    Parameters
    ----------
    rate_constants : dict
        a dictionary of rates.  the keys are tuples of nodes (u,v), the values
        are the rate constants from u to v.
    A, B : iterables
        Groups of nodes specifying the reactant and product groups.  The rates
        returned will be the rate from A to B and vice versa.
    weights : dict
        Dictionary with nodes as keys and weights as values.  The weights are
        the equilibrium occupation probabilities of the nodes in A and B.  They
        are used to do the weighted mean for the final average over inverse
        mean first passage times.
    
    Notes
    -----
    This follows the new graph transformation procedure (NGT) described by 
    David Wales, J. Chem. Phys., 2009 http://dx.doi.org/10.1063/1.3133782
    
    The rate, k_AB computed by this calculation (returned by
    self.get_final_rates) is the inverse mean first passage time averaged over
    the states in A
    
    """
    cdef cNGT* thisptr
    cdef node_list
    cdef node2id
    time_solve = 0.
    def __cinit__(self, rate_constants, A, B, debug=False, weights=None):
        # all this mess is to construct the c++ objects that will be passed to the C++ NGT
        # assign ids to all the nodes
        nodes = set()
        for u, v in rate_constants.iterkeys():
            nodes.add(u)
            nodes.add(v)
        self.node_list = list(nodes)
        self.node2id = dict(( (u, i) for i, u in enumerate(self.node_list) ))
        
        # construct the std::map rate_map
        cdef rate_map_t rate_map
        cdef node_id uid, vid
        for (u, v), k in rate_constants.iteritems():
            uid = self.node2id[u]
            vid = self.node2id[v]
            rate_map[pair_t(uid, vid)] = k
        
        # construct the std::list _A and _B
        cdef stdlist[node_id] _A, _B
        for u in A:
            uid = self.node2id[u]
            _A.push_back(uid)
        for u in B:
            uid = self.node2id[u]
            _B.push_back(uid)
            
        # allocate memory and initialize the NGT object
        self.thisptr = new cNGT(rate_map, _A, _B)
        
        # pass the weights
        cdef map[node_id, double] Peq 
        if weights is not None:
            for u, p in weights.iteritems():
                try:
                    uid = self.node2id[u]
                    Peq[uid] = p
                except KeyError:
                    pass
                        
            self.thisptr.set_node_occupation_probabilities(Peq)
        
        # set the debug flag
        if debug:
            self.thisptr.set_debug()
    
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    
    def compute_rates(self):
        """compute the rates from A->B and B->A"""
        t0 = time.clock()
        self.thisptr.compute_rates()
        self.time_solve = time.clock() - t0 
    
    def compute_rates_and_committors(self):
        """compute the rates from A->B and B->A and the committors for all the intermediates
        
        This is much slower than compute_rates.  Only use this function if you 
        want the committors
        """
        t0 = time.clock()
        self.thisptr.compute_rates_and_committors()
        self.time_solve = time.clock() - t0 
    
    def get_rate_AB(self):
        """return the rate from A->B"""
        return self.thisptr.get_rate_AB()
    
    def get_rate_BA(self):
        """return the rate from B->A"""
        return self.thisptr.get_rate_BA()
    
    def get_rate_AB_SS(self):
        """return the steady state rate from A->B"""
        return self.thisptr.get_rate_AB_SS()
    
    def get_rate_BA_SS(self):
        """return the steady state rate from B->A"""
        return self.thisptr.get_rate_BA_SS()
    
    def get_committors(self):
        """return a dictionary of the committor probabilities"""
        cdef map[node_id, double] qmap = self.thisptr.get_committors() # as reference?
        committors = dict()
        for node, nid in self.node2id.iteritems():
            committors[node] = qmap.at(nid)
        return committors
            
          
    

class NGT(BaseNGT):
    pass
    
