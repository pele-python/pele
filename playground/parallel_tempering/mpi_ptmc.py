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
class MPI_Parallel_Tempering(object):
    """
    Abstract method for MPI Parallel Tempering calculations, it implements all the basic MPI routines. The initialisation
    function and the method to find the swap pattern need to implemented in a specific child method
    """
    __metaclass__  = abc.ABCMeta
    
    def __init__(self, mcrunner, Tmax, Tmin, max_ptiter, pfreq=1):
        self.mcrunner = mcrunner
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size() #total number of processors (replicas)
        self.rank = self.comm.Get_rank() #this is the unique identifier for the process
        self.Tmax = Tmax
        self.Tmin = Tmin
        self.max_ptiter = max_ptiter
        self.ex_outstream = open("exchanges", "w")
        self.ptiter = 0
        self.pfreq = pfreq
        self.no_exchange_int = -12345 #this NEGATIVE number in exchange pattern means that no exchange should be attempted
        self.initialised = False #flag
        self.nodelist = [i for i in xrange(self.nproc)]
        print "processor {0} ready".format(self.rank)
    
    def _scatter_data(self, in_send_array, adim):
        """
        function to scatter data in equal ordered chunks among replica (it relies on the rank of the replica) 
        """
        if (self.rank == 0):
            # process 0 is the root, it has data to scatter
            assert(len(in_send_array) == adim)
            assert(adim % self.nproc == 0) 
            send_array = in_send_array
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
    
    @abc.abstractmethod
    def _find_exchange_buddy(self, Earray):
        """
        This function determines the exchange pattern, this is an abstract methods, it needs to be overwritten.
        An exchange pattern array is constructed, filled with self.no_exchange_int which
        signifies that no exchange should be attempted. This value is replaced with the
        rank of the processor with which to perform the swap if the swap attempt is successful.
        The exchange partner is then scattered to the other processors.
        """
        
    def _exchange_pairs(self, exchange_buddy, data):
        """
        return data from the pair exchange, otherwise return the data unaltered.
        the replica sends to exchange_partner and receives from it (replacing source with self.rank would cause a deadlock)
        """
        if (exchange_buddy != self.no_exchange_int):
            #print "processor {0} p-to-p exchange, old data {1}".format(self.rank, data)
            data = self._point_to_point_exchange_replace(exchange_buddy, exchange_buddy, data) 
            #print "processor {0} p-to-p exchange, new data {1}".format(self.rank, data)
        return data
    
    def _attempt_exchange(self):
        """
        this function brings together all the functions necessary to attempt a configuration swap, it is structures as
        following:
        *root gathers the energies from the slaves
        
        """
        #gather energies, only root will do so
        Earray = self._gather_energies(self.energy)
        if Earray is not None:
            print "Earray", Earray
        #find exchange buddy
        exchange_buddy = self._find_exchange_buddy(Earray)
        #attempt configurations swap
        self.config = self._exchange_pairs(exchange_buddy, self.config)
        #swap energies
        E = self._exchange_pairs(exchange_buddy, np.array([self.energy],dtype='d'))
        assert(len(E)==1)
        self.energy = E[0]
    
    @abc.abstractmethod
    def _initialise(self):
        """
        perform all the tasks required prior to starting the computation
        """
    
    @abc.abstractmethod
    def _print(self):
        """this function is responsible for printing and/or dumping the data, let it be printing the histograms or else"""
    
    def one_iteration(self):
        """Perform one parallel tempering iteration, this consists of the following steps:
        *set the coordinates
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
            print "processor {0} iteration {1}".format(self.rank,ptiter)
            self.one_iteration()
            ptiter += 1

class MPI_PT_Simple(MPI_Parallel_Tempering):
    """
    This class performs parallel tempering by alternating swaps with right and left neighbours with geometrically 
    distributed temperatures.
    """
    def __init__(self, mcrunner, Tmax, Tmin, max_ptiter=10, pfreq=1):
        super(MPI_PT_Simple,self).__init__(mcrunner, Tmax, Tmin, max_ptiter, pfreq=pfreq)
        self.exchange_dic = {1:'right',-1:'left'}
        self.exchange_choice = random.choice(self.exchange_dic.keys()) 
        self.anyswap = False #set to true if any swap will happen
        self.permutation_pattern = np.zeros(self.nproc,dtype='int32') #this is useful to print exchange permutations
        
    def _print(self):
        base_directory = "ptmc_results"
        if (self.ptiter == 0 and self.rank == 0):
            if not os.path.exists(base_directory):
                os.makedirs(base_directory)
            fname = "{0}/temperatures".format(base_directory)
            np.savetxt(fname, self.Tarray, delimiter='\t')
        
        if (self.rank == 0 and self.anyswap == True):
            fname = "{0}/rem_permutations".format(base_directory)
            f = open(fname,'a')
            iteration = self.mcrunner.get_iterations_count()
            f.write('{0}\t'.format(iteration))
            for p in self.permutation_pattern:
                f.write('{0}\t'.format(p))
            f.write('\n')
            f.close()
                        
        if (self.ptiter % self.pfreq == 0):
            directory = "{0}/{1}".format(base_directory,self.rank)
            if not os.path.exists(directory):
                os.makedirs(directory)
            iteration = self.mcrunner.get_iterations_count()
            fname = "{0}/Visits.his.{1}".format(directory,float(iteration))
            self.mcrunner.dump_histogram(fname)
        
        if (self.ptiter == self.max_ptiter-1):
            directory = "{0}/{1}".format(base_directory,self.rank)
            fname = "{0}/parameters".format(directory)
            accepted_frac = self.mcrunner.get_accepted_fraction()
            init_stepsize = self.mcrunner.stepsize
            fin_stepsize = self.mcrunner.get_stepsize()
            ncount = self.mcrunner.get_iterations_count()
            f = open(fname,'a')
            f.write('node:\t{0}\n'.format(self.rank))
            f.write('temperature:\t{0}\n'.format(self.T))
            f.write('initial step size:\t{0}\n'.format(init_stepsize))
            f.write('adapted step size:\t{0}\n'.format(fin_stepsize))
            f.write('PT iterations:\t{0}\n'.format(self.ptiter))
            f.write('total MC iterations:\t{0}\n'.format(ncount))
            f.write('acceptance fraction:\t{0}\n'.format(accepted_frac))
            f.close()
            
    
    def _get_temps(self):
        """
        set up the temperatures by distributing them exponentially. We give root the lowest temperature.
        This should increase performance when pair lists are used (they are updated less often at low temperature
        or when steps involve minimisation, as the low temperatures are closer to the minimum)
        """
        if (self.rank == 0):
            CTE = np.exp( np.log( self.Tmax / self.Tmin ) / (self.nproc-1) )
            Tarray = [self.Tmin * CTE**i for i in range(self.nproc)]
            #Tarray = np.linspace(self.Tmin,self.Tmax,num=self.nproc)
            self.Tarray = np.array(Tarray[::-1],dtype='d')
        else:
            self.Tarray = None
    
    def _initialise(self):
        """
        perform all the tasks required prior to starting the computation
        """
        self._get_temps()
        self.config, self.energy = self.mcrunner.get_config()
        self.T = self._scatter_single_value(self.Tarray)
        print "processor {0} temperature {1}".format(self.rank,self.T)
        self.mcrunner.set_control(self.T)
        self.initialised = True
    
    def _find_exchange_buddy(self, Earray):
        """
        This function determines the exchange pattern alternating swaps with right and left neighbours.
        An exchange pattern array is constructed, filled with self.no_exchange_int which
        signifies that no exchange should be attempted. This value is replaced with the
        rank of the processor with which to perform the swap if the swap attempt is successful.
        The exchange partner is then scattered to the other processors.
        """        
        if (self.rank == 0):
            assert(len(Earray)==len(self.Tarray))
            exchange_pattern = np.empty(len(Earray),dtype='int32')
            exchange_pattern.fill(self.no_exchange_int)
            self.anyswap = False
            for i in xrange(0,self.nproc,2):
                print 'exchange choice: ',self.exchange_dic[self.exchange_choice] #this is a print statement that has to be removed after initial implementation
                E1 = Earray[i]
                T1 = self.Tarray[i]
                E2 = Earray[i+self.exchange_choice]
                T2 = self.Tarray[i+self.exchange_choice]
                deltaE = E1 - E2
                deltabeta = 1./T1 - 1./T2
                w = min( 1. , np.exp( deltaE * deltabeta ) )
                rand = np.random.rand()
                #print "E1 {0} T1 {1} E2 {2} T2 {3} w {4}".format(E1,T1,E2,T2,w) 
                if w > rand:
                    #accept exchange
                    self.ex_outstream.write("accepting exchange %d %d %g %g %g %g %d\n" % (self.nodelist[i], self.nodelist[i+self.exchange_choice], E1, E2, T1, T2, self.ptiter))
                    assert(exchange_pattern[i] == self.no_exchange_int)                      #verify that is not using the same processor twice for swaps
                    assert(exchange_pattern[i+self.exchange_choice] == self.no_exchange_int) #verify that is not using the same processor twice for swaps
                    exchange_pattern[i] = self.nodelist[i+self.exchange_choice]
                    exchange_pattern[i+self.exchange_choice] = self.nodelist[i]
                    self.anyswap = True
            ############end of for loop###############
            #record self.permutation_pattern to print permutations in print function
            if self.anyswap:
                for i,buddy in enumerate(exchange_pattern):
                    if (buddy != self.no_exchange_int):
                        self.permutation_pattern[i] = buddy+1 #to conform to fortran notation
                    else:
                        self.permutation_pattern[i] = i+1 #to conform to fortran notation
        else:
            exchange_pattern = None
        
        self.exchange_choice *= -1 #swap direction of exchange choice
        #print "exchange_pattern",exchange_pattern
        #now scatter the exchange pattern so that everybody knows who their buddy is
        exchange_buddy = self._scatter_single_value(np.array(exchange_pattern,dtype='d'))
        #print "processor {0} buddy {1}".format(self.rank,int(exchange_buddy))
        return int(exchange_buddy)
            
            
            
    
    
    
        
                
            
              
                
                
                