import numpy as np
import multiprocessing as mp

import pygmin.utils.fix_multiprocessing

from pygmin.landscape import DoubleEndedConnect
from pygmin.landscape.connect_min import _refineTS
from pygmin.transition_states import NEBPar

__all__ = ["DoubleEndedConnectPar"]

def _refineTSWrapper(inputs):
    """
    a stupid wrapper to allow _refineTS to be used with Pool.map which only supports one argument
    """
    return _refineTS(inputs[0], inputs[1], **inputs[2])


class DoubleEndedConnectPar(DoubleEndedConnect):
    """
    Overload some of the routines from DoubleEndedConnect so they can
    be parallelized
    
    Extra Parameters
    ----------------
    ncores :
        the number of cores to use in parallel runs
        
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
        return super(DoubleEndedConnectPar, self).__init__(*args, **kwargs)

    def _refineTransiitonStates(self, neb, climbing_images):
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
            print ""
            print "refining transition state from NEB climbing image:", count, "out of", nrefine
            coords = neb.coords[i,:]
            #get guess for initial eigenvector from NEB tangent
            if True:
                eigenvec0 = neb.tangent( [neb.energies[i-1], neb.coords[i-1,:]],
                                         [neb.energies[i], neb.coords[i,:]],
                                         [neb.energies[i+1], neb.coords[i+1,:]]
                                        )
            input_args.append((self.pot, coords, {"tsSearchParams":self.tsSearchParams, "eigenvec0":eigenvec0}))

        #do all the transition state searches in parallel using a pool of workers
        print "refining transition states in parallel on", self.ncores, "cores"
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
                    goodts = self._addTransitionState(tsret.energy, tsret.coords, m1ret, m2ret, tsret.eigenvec, tsret.eigenval)
                    if goodts:
                        ngood_ts += 1
        except:
            #It's important to make sure the child processes are closed even
            #if when an exception is raised.  
            #Note: I'm not sure this is the right way to do it.
            print "exception raised while running multiprocessing.Pool.map"
            print "  terminating pool"
            mypool.terminate()
            mypool.join()
            raise
        mypool.close()
        mypool.join()
        
        print "found", ngood_ts, "good transition states from", nrefine, "candidates"
        return ngood_ts > 0
    
    def _getNEB(self, *args, **kwargs):
        return NEBPar(*args, ncores=self.ncores, **kwargs)



if __name__ == "__main__":
    from pygmin.landscape.connect_min import test
    test(DoubleEndedConnectPar)
