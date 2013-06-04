import multiprocessing as mp
import logging

#this import fixes some bugs in how multiprocessing deals with exceptions
import pele.utils.fix_multiprocessing

from pele.landscape import DoubleEndedConnect, LocalConnect
from pele.landscape.local_connect import _refineTS
from pele.transition_states import create_NEB

__all__ = ["DoubleEndedConnectPar", "LocalConnectPar"]

logger = logging.getLogger("pele.connect")

def _refineTSWrapper(inputs):
    """
    a stupid wrapper to allow _refineTS to be used with Pool.map which only supports one argument
    """
    return _refineTS(inputs[0], inputs[1], **inputs[2])


class DoubleEndedConnectPar(DoubleEndedConnect):
    """
    parallelized version of DoubleEndedConnect
    
    Overload some of the routines from DoubleEndedConnect so they can
    be parallelized
    
    Parameters
    ----------------
    ncores : int, optional
        the number of cores to use in parallel runs
        
    Notes
    -----
    This class inherits from DoubleEndedConnect, so can accepts all those parameters as well.
    
    The routines that are done in parallel are::
    
    1. NEB : the potentials for each image are calculated in parallel
    2. findTransitionStates : each transition state candidate from the NEB run is refined in parallel. 
    
    See Also
    --------
    DoubleEndedConnect : the class this inherits from
    LocalConnectPar : parallel version of LocalConnect
    
    """
    def __init__(self, *args, **kwargs):
        #self.ncores = ncores
        try:
            self.ncores = kwargs.pop("ncores")
        except KeyError:
            self.ncores = 4
        return super(DoubleEndedConnectPar, self).__init__(*args, **kwargs)

    def _getLocalConnectObject(self):
        return LocalConnectPar(self.pot, self.mindist, ncores=self.ncores, **self.local_connect_params)


class LocalConnectPar(LocalConnect):
    """
    Overload some of the routines from LocalConnect so they can
    be parallelized
    
    Parameters
    ----------------
    inherited params :
        all required and optional parameters from LocalConnect are also accepted
    ncores :
        the number of cores to use in parallel runs
    
    See Also
    --------
    LocalConnect : base class
    DoubleEndedConnectPar : uses this class
    pele.transition_states.NEBPar : parallel version of NEB

    
    Notes
    -----    
    The routines that are done in parallel are:
    
    NEB : the potentials are calculated in parallel
    
    findTransitionStates : each transition state candidate from the NEB
        run is refined in parallel. 
    
    """
    def __init__(self, *args, **kwargs):
        #self.ncores = ncores
        try:
            self.ncores = kwargs.pop("ncores")
        except KeyError:
            self.ncores = 4
        super(LocalConnectPar, self).__init__(*args, **kwargs)
        self.NEBparams["parallel"] = True
        self.NEBparams["ncores"] = self.ncores

    def _refineTransitionStates(self, neb, climbing_images):
        """
        Use a pool of parallel workers to calculate the transition states.
        """
        #find the nearest transition state to the transition state candidates
        nrefine = min(self.nrefine_max, len(climbing_images))
        if nrefine == 0:
            return False
        count = 0
        #prepare the input arguments for all the transition state searches
        input_args = []
        for energy, i in climbing_images[:nrefine]:
            count += 1
            logger.info("\nrefining transition state from NEB climbing image: %s %s %s", count, "out of", nrefine)
            coords = neb.coords[i,:]
            #get guess for initial eigenvector from NEB tangent
            if True:
                eigenvec0 = neb.tangent( neb.energies[i], neb.energies[i-1], neb.energies[i+1],
                                         neb.distance(neb.coords[i,:], neb.coords[i-1,:])[1],
                                         neb.distance(neb.coords[i,:], neb.coords[i+1,:])[1],
                                        )

            input_args.append((self.pot, coords, {"tsSearchParams":self.tsSearchParams, "eigenvec0":eigenvec0}))

        #do all the transition state searches in parallel using a pool of workers
        logger.info("refining transition states in parallel on %s %s", self.ncores, "cores")
        mypool = mp.Pool(self.ncores)
        try:
            #there is a bug in Python so that exceptions in multiprocessing.Pool aren't
            #handled correctly.  A fix is to add a timeout (.get(timeout))
            returnlist = mypool.imap_unordered( _refineTSWrapper, input_args )#.get(9999999)
            
            ngood_ts = 0
            #analyze the return values and 
            #find the minimum on either side of each good transition state
            for ret in returnlist:#.get(99999999):        
                ts_success = ret[0]
                if ts_success:
                    #the transition state is good, add it to the graph
                    tsret, m1ret, m2ret = ret[1:4]
                    self.res.new_transition_states.append( (tsret, m1ret, m2ret) )
                    ngood_ts += 1
        except:
            #It's important to make sure the child processes are closed even
            #if when an exception is raised.  
            #Note: I'm not sure this is the right way to do it.
            logger.error("exception raised while running multiprocessing.Pool.map")
            logger.error("  terminating pool")
            mypool.terminate()
            mypool.join()
            raise
        mypool.close()
        mypool.join()
        
        logger.info("found %s %s %s %s", ngood_ts, "good transition states from", nrefine, "candidates")
        return ngood_ts > 0
    
#    def _getNEB(self, *args, **kwargs):
#        #this is all that need be changed to get the NEB to run in parallel.
#        return create_NEB(*args, parallel=True, ncores=self.ncores, **kwargs)
##        return NEBPar(*args, ncores=self.ncores, **kwargs)



if __name__ == "__main__":
    from pele.landscape.connect_min import test
    test(DoubleEndedConnectPar, natoms=28)
