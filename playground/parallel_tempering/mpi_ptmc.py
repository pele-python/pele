from __future__ import division
import numpy as np
from mpi4py import MPI

class MPI_Parallel_Tempering(object):
    def __init__(self, mc_runner, Tmax, Tmin, mciter = 100):
        self.mcrunner = mc_runner
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size() #total number of processors (replicas)
        self.rank = self.comm.Get_rank() #this is the unique identifier for the process  
        self.Tmax = Tmax
        self.Tmin = Tmin    
        
    def _get_temps(self):
        """
        set up the temperatures
        distribute them exponentially
        """
        CTE = np.exp( np.log( self.Tmax / self.Tmin ) / (self.nproc-1) )
        Tarray = np.array([self.Tmin * CTE**i for i in range(self.nproc)],dtype='double')
        return Tarray
    
    def _scatter_data(self, in_send_array):
        """
        function to scatter data in equal ordered chunks among replica (it relies on the rank of the replica) 
        """
        assert(len(in_send_array) % self.nproc == 0) 
        
        if (self.rank == 0):
            send_array = in_send_array
        else:
            send_array = None
        
        recv_array = np.zeros(len(send_array) /self.nproc,dtype='double')
        self.comm.Scatter(send_array, recv_array, root=0) 
        return recv_array 
    
    def _scatter_temps(self, Tarray):
        """
        returns the assigned temperature for each replica
        """
        T = self._scatter_data(Tarray)
        assert(len(T) == 1)
        return T[0]
    
    def _broadcast_data(self, in_data):
        """
        identical data are broadcasted from root to all other processes
        """
        if(self.rank == 0):
            bcast_data = in_data
        else:
            bcast_data = None
        bcast_data = self.comm.Bcast(bcast_data, root=0)
        return bcast_data
    
    def _gather_data(self, in_send_array):
        """
        function to gather data in equal ordered chunks from replicas (it relies on the rank of the replica) 
        """
        if (self.rank == 0):
            recv_array = np.zeros(len(in_send_array) * self.nproc,dtype='double')
        else:
            recv_array = None
        
        self.comm.Gather(in_send_array, recv_array, root=0)
        
        if (self.rank != 0):
            assert(recv_array == None)
        
        return recv_array
    
    def _gather_energies(self, E):
        send_Earray = np.array([E],dtype='double')
        recv_Earray = self._gather_data(send_Earray)
        return recv_Earray
    
    def _point_to_point_exchange_replace(self, dest, source, data):
        assert(self.rank == dest or self.rank == source)
        self.comm.Sendrecv_replace(data, dest=dest,source=source)
        return data
        
    
    def tryExchangePar(self):
        #choose which pair to try and exchange
        k = np.random.random_integers( 0, self.nreplicas - 2)
        #print "trying exchange", k, k+1
        #determine if the exchange will be accepted
        E1, T1 = self.getRepEnergyT(k)
        E2, T2 = self.getRepEnergyT(k+1)
        deltaE = E1 - E2
        deltabeta = 1./T1 - 1./T2
        w = min( 1. , np.exp( deltaE * deltabeta ) )
        rand = np.random.rand()
        if w > rand:
            #accept exchange
            self.ex_outstream.write("accepting exchange %d %d %g %g %g %g %d\n" % (k, k+1, E1, E2, T1, T2, self.step_num) )
            self.doExchangePar(k, k+1)
        #else:
            #self.ex_outstream.write("rejecting exchange %d %d %g %g\n" % (k, k+1, w, rand) )
            #print "rejecting exchange ", k, k+1, w, rand