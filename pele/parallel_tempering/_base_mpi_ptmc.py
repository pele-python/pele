from __future__ import division
import abc
import numpy as np
import random
import os
from mpi4py import MPI

"""
An optimal Parallel Tempering strategy should make sure that all MCMC walks take roughly the same amount of time. 
Besides this fundamental consideration, note that root (rank=0) is not an evil master but rather an enlightened dictator that leads by
example: root is responsible to assign jobs and control parameters (e.g. temperature) to the slaves but it also performs MCMC walks along 
with them. For this reason it might be optimal to give root a set of control parameters for which the simulation is leaner so that it 
can start doing its own things while the slaves finish their work.  
"""
class _MPI_Parallel_Tempering(object):
    """
    Abstract method for MPI Parallel Tempering calculations, it implements all the basic MPI routines. The initialisation
    function and the method to find the swap pattern need to implemented in a specific child method
    """
    __metaclass__  = abc.ABCMeta
    
    def __init__(self, mcrunner, Tmax, Tmin, max_ptiter, pfreq=1, base_directory=None, verbose=False):
        self.mcrunner = mcrunner
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size() #total number of processors (replicas)
        self.rank = self.comm.Get_rank() #this is the unique identifier for the process
        self.Tmax = Tmax
        self.Tmin = Tmin
        self.max_ptiter = max_ptiter
        self.ex_outstream = open("exchanges", "w")
        self.verbose = verbose
        self.ptiter = 0
        self.pfreq = pfreq
        self.no_exchange_int = -12345 #this NEGATIVE number in exchange pattern means that no exchange should be attempted
        self.initialised = False #flag
        self.nodelist = [i for i in xrange(self.nproc)]
        self.swap_accepted_count = 0
        self.swap_rejected_count = 0
        if base_directory is None:
            self.base_directory = os.path.join(os.getcwd(),'ptmc_results')
        else:
            self.base_directory = base_directory
        
        print "processor {0} ready".format(self.rank)
    
    @abc.abstractmethod
    def _find_exchange_buddy(self, Earray):
        """
        This function determines the exchange pattern, this is an abstract methods, it needs to be overwritten.
        An exchange pattern array is constructed, filled with self.no_exchange_int which
        signifies that no exchange should be attempted. This value is replaced with the
        rank of the processor with which to perform the swap if the swap attempt is successful.
        The exchange partner is then scattered to the other processors.
        """
    
    @abc.abstractmethod
    def _initialise(self):
        """
        perform all the tasks required prior to starting the computation, including initialising the output files
        """
    
    @abc.abstractmethod
    def _print(self):
        """this function is responsible for printing and/or dumping the data, let it be printing the histograms or else"""
    
    def one_iteration(self):
        """Perform one parallel tempering iteration, this consists of the following steps:
        *set the coordinates#now scatter the exchange pattern so that everybody knows who their buddy is
        *run the MCrunner for a predefined number of steps
        *collect the results (energy and new coordinates)
        *attempt an exchange
        """
        #set configuration and temperature at which want to perform run
        self.mcrunner.set_config(self.config, self.energy)
        #now run the MCMC walk
        self.mcrunner.run()
        #collect the results
        result = self.mcrunner.get_results()
        self.energy = result.energy
        self.config = result.coords
        self._attempt_exchange()
        #print and increase parallel tempering count
        self._print()
        self.ptiter += 1
         
            
    def run(self):
        """Run consists of multiple single iterations, plus initialisation if MPI_PT has not been initialised yet 
        """
        if self.initialised is False:
            self._initialise()
        ptiter = 0
        while (ptiter < self.max_ptiter):
            if self.verbose:
                print "processor {0} iteration {1}".format(self.rank,ptiter)
            self.one_iteration()
            ptiter += 1
    
    def _scatter_data(self, in_send_array, adim):
        """
        function to scatter data in equal ordered chunks among replica (it relies on the rank of the replica)
        note that arrays of doubles are scattered 
        """
        if (self.rank == 0):
            # process 0 is the root, it has data to scatter
            assert(len(in_send_array) == adim)
            assert(adim % self.nproc == 0) 
            send_array = np.array(in_send_array,dtype='d')
        else:
            # processes other than root do not send
            assert(adim % self.nproc == 0) 
            send_array = None
        
        recv_array = np.empty(adim/self.nproc,dtype='d')
        self.comm.Scatter(send_array, recv_array, root=0) 
        return recv_array 
    
    def _scatter_single_value(self, send_array):
        """
        returns a single value from a scattered array for each replica (e.g. Temperature or swap partner)
        this implies that send array must be of the same length as the number of processors
        """
        if (self.rank == 0):
            assert(len(send_array) == self.nproc)
        
        T = self._scatter_data(send_array, self.nproc)
        assert(len(T) == 1)
        return T[0]
    
    def _broadcast_data(self, in_data, adim):
        """
        identical data are broadcasted from root to all other processes
        """
        if(self.rank == 0):
            bcast_data = in_data
        else:
            bcast_data = np.empty(adim,dtype='d')
        bcast_data = self.comm.Bcast(bcast_data, root=0)
        return bcast_data
    
    def _gather_data(self, in_send_array):
        """
        function to gather data in equal ordered chunks from replicas (it relies on the rank of the replica)
        note that gather assumes that all the subprocess are sending the same amount of data to root, to send
        variable amounts of data must use the MPI_gatherv directive 
        """
        if (self.rank == 0):
            recv_array = np.zeros(len(in_send_array) * self.nproc,dtype='d')
        else:
            recv_array = None
        
        self.comm.Gather(in_send_array, recv_array, root=0)
        
        if (self.rank != 0):
            assert(recv_array == None)
        
        return recv_array
    
    def _gather_energies(self, E):
        """
        gather energy configurations from all processors
        """
        send_Earray = np.array([E],dtype='d')
        recv_Earray = self._gather_data(send_Earray)
        return recv_Earray
    
    def _point_to_point_exchange_replace(self, dest, source, data):
        """
        swap data between two processors, the message sent buffer is replaced with the received message
        """
        assert(dest == source)
        self.comm.Sendrecv_replace(data, dest=dest,source=source)
        return data
    
    def _exchange_pairs(self, exchange_buddy, data):
        """
        return data from the pair exchange, otherwise return the data unaltered.
        the replica sends to exchange_partner and receives from it (replacing source with self.rank would cause a deadlock)
        """
        if (exchange_buddy != self.no_exchange_int):
            #print "processor {0} p-to-p exchange, old data {1}".format(self.rank, data)
            data = self._point_to_point_exchange_replace(exchange_buddy, exchange_buddy, data) 
            #print "processor {0} p-to-p exchange, new data {1}".format(self.rank, data)
            self.swap_accepted_count+=1
        else:
            self.swap_rejected_count+=1
        return data
    
    def _attempt_exchange(self):
        """
        this function brings together all the functions necessary to attempt a configuration swap, it is structures as
        following:
        *root gathers the energies from the slaves
        
        """
        #gather energies, only root will do so
        Earray = self._gather_energies(self.energy)
        if self.verbose:
            if Earray is not None:
                print "Earray", Earray
        #find exchange pattern (list of exchange buddies)
        exchange_pattern = self._find_exchange_buddy(Earray)
        #now scatter the exchange pattern so that everybody knows who their buddy is
        exchange_buddy = self._scatter_single_value(np.array(exchange_pattern,dtype='d'))
        exchange_buddy = int(exchange_buddy)
        #attempt configurations swap
        self.config = self._exchange_pairs(exchange_buddy, self.config)
        #swap energies
        E = self._exchange_pairs(exchange_buddy, np.array([self.energy],dtype='d'))
        assert(len(E)==1)
        self.energy = E[0]
            
            
            
    
    
    
        
                
            
              
                
                
                