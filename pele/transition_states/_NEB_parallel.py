import numpy as np
import copy
import multiprocessing as mp
import logging

from pele.transition_states import NEB

__all__ = ["NEBPar"]

logger = logging.getLogger("pele.connect.neb")

class _PotentialProcess(mp.Process):
    """
    this class defines the worker object to be run as a separate process
    
    Parameters
    ----------
    pot : 
        the potential object
    conn : 
        the communications pipe
    copy_potential : integer
        if copy_potential > 1 then the copy_potential deep copies of the potential
        will be made and each one used for a different item in coordslist.  This
        is useful to minimize rebuilding of neighborlists.

    Notes
    -----
    attributes:
    
    start():
        start the process.  The proceess will sit idle waiting for instructions
        sent through the pipe conn
    join():
        wait until the process has finished it's task.  (Send message "kill" to stop 
        the process waiting for the next task)
    terminate():
        kill the process
    
    
    Examples
    ---------
    example usage::
    
    >>> parent_conn, child_conn = mp.Pipe()
    >>> worker = _PotentialProcess( potential, child_conn)
    >>> worker.start()
    
    >>> X1 = np.random.rand(3*10)
    >>> X2 = np.random.rand(3*10)
    >>> parent_conn.send(("calculate energy gradient", [X1, X2]))
    >>> eglist = parent_conn.recv()
    >>> print "energies", eglist[0][0], eglist[1][0]
    
    >>> parent_conn.send(("kill",))
    >>> worker.join()
    >>> worker.terminate() 


    
    once this process is started (start()) it will be waiting for a message for what to do
    """
    def __init__(self, pot, conn, copy_potential=0):
        mp.Process.__init__(self)
        self.conn = conn
        self.copy_potential = int(copy_potential)
        if self.copy_potential > 1:
            self.potlist = [copy.deepcopy(pot) for i in range(self.copy_potential)] 
        else:
            self.pot = pot
    
    def getEnergyGradientMultiple(self, coordslist):
        eglist = []
        if self.copy_potential > 1:
            assert len(coordslist) == self.copy_potential
            for i, coords in enumerate(coordslist):
                eglist.append( self.potlist[i].getEnergyGradient(coords) )
        else:
            for coords in coordslist:
                eglist.append( self.pot.getEnergyGradient(coords) )
        return eglist

    def run(self):
        #this redefines mp.Process.run
        while 1:
            message = self.conn.recv()
            #print "message", message
            if message[0] == "kill":
                #print "terminating", self.name
                return
            elif message[0] == "calculate energy gradient":
                #print "received message: calculating gradient"
                coordslist = message[1]
                eglist = self.getEnergyGradientMultiple(coordslist)
                self.conn.send(eglist)
            else:
                logger.error("unknown message: %s\n%s", self.name, message)
                
                
class NEBPar(NEB):
    """
    Exactly the same as NEB, but calculate the potentials in parallel
    
    Parameters
    --------------
    inherited parameters :
        all required and optional parameters from NEB are accepted
    ncores : 
        the number of cores to use
    
    See Also
    --------
    NEB : base class
    pele.landscape.LocalConnectPar : were this class is used
    """
    def __init__(self, *args, **kwargs):
        #self.ncores = ncores
        try:
            self.ncores = kwargs.pop("ncores")
        except KeyError:
            self.ncores = 4
        #we must deal with copy_potential separately, so extract copy potential
        #from the parameter list and set it to false.
        try:
            self.par_copy_potential = kwargs.pop("copy_potential")
        except KeyError:
            self.par_copy_potential = False
        kwargs["copy_potential"] = False 
            
        ret = super(NEBPar, self).__init__(*args, **kwargs)

        #we need to distribute self.nimages potential calculations over self.ncores processors
        #if there are 5 jobs and 3 processors it will be distributed like this
        #core1 will have images 1,2 
        #core2 will have images 3,4 
        #core3 will have image  5 
        self.par_images = [] #a list of lists of which images are on which cores
        nall = int(self.nimages / self.ncores) #all cores have at least this many images
        nextra = self.nimages - nall * self.ncores
        count = 0
        for i in range(self.ncores):
            if i < nextra:
                #this core has nall+1 images
                self.par_images.append( np.array(range(count,count+nall+1)))
                count += nall + 1
            else:
                #this core has nall images
                self.par_images.append( np.array(range(count,count+nall)))
                count += nall
        #print self.par_images        
        
        #initialize the parallel workers with the potential
        self.workerlist = []
        self.connlist = []
        for i in range(self.ncores):
            parent_conn, child_conn = mp.Pipe()
            worker_nimages = len(self.par_images[i])
            if not self.par_copy_potential:
                worker_nimages = 0
            worker = _PotentialProcess( self.potential, child_conn, copy_potential=worker_nimages)
            self.workerlist.append( worker )
            self.connlist.append( parent_conn )
            #worker.start()

        return ret

    def _getRealEnergyGradient(self, coordsall):
        """
        override the function from NEB to calculate the energies in parallel
        """
        #start all workers 
        for i in range(self.ncores):
            coordslist = coordsall[self.par_images[i],:]
            self.connlist[i].send(("calculate energy gradient", coordslist))
        
        #receive energies and gradients from all workers
        realgrad = np.zeros(coordsall.shape)
        for i in range(self.ncores):
            eglist = self.connlist[i].recv()
            imagelist = self.par_images[i]
            for j in range(len(imagelist)):
                image = imagelist[j]
                self.energies[image] = eglist[j][0]
                realgrad[image,:]  = eglist[j][1]
            
        return realgrad   
    
    def _killWorkers(self):
        for i, worker in enumerate(self.workerlist):
            #print "killing"
            #tell the worker to stop waiting for messages
            self.connlist[i].send(("kill",))
            #print "joining"
            #wait for the worker to finish what it's doing
            worker.join()
            #kill the worker process
            worker.terminate()
            worker.join()
            #print "is alive?", worker.is_alive()


    def optimize(self, *args, **kwargs):
        """
        wrap the optimize routine of NEB so the workers can be started
        and, importantly, stopped, even if an exception is raised.
        """
        try:
            for worker in self.workerlist:
                worker.start()
            logger.info("running NEB in parallel with %s %s", self.ncores, "cores")
            ret = super(NEBPar, self).optimize(*args, **kwargs)
            self._killWorkers()
            return ret
        except:
            logger.error("exception raised while doing NEB in parallel, terminating child processes")
            for worker in self.workerlist:
                worker.terminate()
                worker.join()
            #self._killWorkers()
            raise
    
if __name__ == "__main__":
    logger.basicConfig(level=logger.DEBUG)
    from pele.transition_states._NEB import nebtest
    nebtest(NEBPar)

    
