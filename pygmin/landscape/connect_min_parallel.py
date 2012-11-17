import numpy as np
import multiprocessing as mp

from pygmin.landscape import DoubleEndedConnect
from pygmin.landscape.connect_min import _refineTS
from pygmin.transition_states import NEBPar

__all__ = ["DoubleEndedConnectPar"]

def _refineTSWrapper(inputs):
    """
    a stupid wrapper to allow _refineTS to be used with map which only allows one argument
    """
    return _refineTS(inputs[0], inputs[1], **inputs[2])


class DoubleEndedConnectPar(DoubleEndedConnect):
    def __init__(self, *args, **kwargs):
        #self.ncores = ncores
        try:
            self.ncores = kwargs.pop("ncores")
        except KeyError:
            self.ncores = 4
        return super(DoubleEndedConnectPar, self).__init__(*args, **kwargs)

    def _refineTransiitonStates(self, neb, climbing_images):
        """
        start each transition state search in a separate process, up to a max of
        NPAR.
        """
        #find the nearest transition state to the transition state candidates
        nrefine = min(self.nrefine_max, len(climbing_images))
        if nrefine == 0:
            return False
        count = 0
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

        print "refining transition states in parallel on", self.ncores, "cores"
        mypool = mp.Pool(self.ncores)
        try:
            returnlist = mypool.map( _refineTSWrapper, input_args )
        except:
            #It's important to make sure the child processes are closed even
            #if when an exception is raised.  
            #Note: I'm not sure this is the right way to do it.
            print "exception raised while running multiprocessing.Pool.map"
            print "  closing pool"
            mypool.close()
            print "  joining pool"
            mypool.join()
            raise
        mypool.close()
        mypool.join()
        ngood_ts = 0
        for ret in returnlist:        
            ts_success = ret[0]
            if ts_success:
                #the transition state is good, add it to the graph
                tsret, m1ret, m2ret = ret[1:4]
                goodts = self._addTransitionState(tsret.energy, tsret.coords, m1ret, m2ret, tsret.eigenvec, tsret.eigenval)
                if goodts:
                    ngood_ts += 1
        print "found", ngood_ts, "good transition states from", nrefine, "candidates"
        return ngood_ts > 0
    
    def _getNEB(self, *args, **kwargs):
        """
        wrap the actual call to NEB so it can be overloaded
        """
        return NEBPar(*args, ncores=self.ncores, **kwargs)



if __name__ == "__main__":
    from pygmin.landscape.connect_min import test
    test(DoubleEndedConnectPar)
