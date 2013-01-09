from pygmin.optimize import Result
import pygmin.defaults as defaults
from pygmin.transition_states import NEB, InterpolatedPath, findTransitionState, minima_from_ts, create_NEB

__all__ = ["LocalConnect"]

def _refineTS(pot, coords, tsSearchParams=dict(), eigenvec0=None):
    """
    find nearest transition state to NEB climbing image.  Then fall
    off the transition state to find the associated minima.
    
    This would naturally be a part of DoubleEndedConnect.  I separated it
    to make it more easily parallelizable.      
    """
    #run ts search algorithm
    kwargs = dict(defaults.tsSearchParams.items() + tsSearchParams.items())
    ret = findTransitionState(coords, pot, eigenvec0=eigenvec0, **kwargs)
    
    #check to make sure it is a valid transition state 
    coords = ret.coords
    if not ret.success:
        print "transition state search failed"
        return False, None
        
    if ret.eigenval >= 0.:
        print "warning: transition state has positive lowest eigenvalue, skipping:", ret.eigenval, ret.energy, ret.rms
        print "         not adding transition state"
        return False, None
    
    #find the minima which this transition state connects
    print "falling off either side of transition state to find new minima"
    ret1, ret2 = minima_from_ts(pot.getEnergyGradient, coords, n = ret.eigenvec, \
        displace=1e-3, quenchParameters={"tol":1e-7, "iprint":-1})
    
    return True, ret, ret1, ret2



class LocalConnect(object):
    """
    a class to do a single local connect run, i.e. NEB + transition state search

    Parameters
    ----------
    pot : potential object
        the potential
    mindist : callable
        the function which returns the optimized minimum distance between
        two structures
    tsSearchParams: dict
        parameters passed to the transition state search algorithm
    NEBparams : dict
        NEB setup parameters.  Use NEBquenchParams for parameters related 
        to the optimization of the band.
    NEBquenchParams : dict
        parameters passed to the NEB minimization routine
    nrefine_max : int
        the maximum number of NEB transition state candidates to refine
    reoptimize_climbing : int
        the number of iterations to use for re-optimizing the climbing images
        after the NEB is done.
    verbosity : int
        this controls how many status messages are printed.  (not really
        implemented yet)
    
    
    Notes
    -----
    this class takes two minima as input, does an NEB run to find transition state
    candidates, then refines those candidates into transition states.
    
    This is the core routine of DoubleEndedConnect.  It is separated out in order
    to make parallelization easier.  This class intentionally has no knowledge of the
    global landscape (database, graph, etc.).
    
    See Also
    --------
    DoubleEndedConnect : the routine from which local connect is genearlly called
    pygmin.transition_states.NEB : one of the core routines
    pygmin.transition_states.create_NEB : the wrapper which sets up NEB
    pygmin.transition_states.findTransitionState : one of the core routine
    LocalConnectPar : parallel version of this class
    
    """
    def __init__(self, pot, mindist, tsSearchParams=dict(), 
                 NEBquenchParams = dict(), verbosity=1,
                 NEB_image_density = 10., NEB_iter_density=15., NEBparams=dict(), 
                 nrefine_max=100, reoptimize_climbing=0, 
                 NEB_max_images=40):
        self.pot = pot
        self.mindist = mindist
        self.tsSearchParams = tsSearchParams
        self.NEBquenchParams = NEBquenchParams
        self.verbosity = int(verbosity)
        self.nrefine_max = nrefine_max
        
        self.NEB_image_density = float(NEB_image_density)
        self.NEB_iter_density = float(NEB_iter_density)
        self.NEBparams = NEBparams
        self.reoptimize_climbing = reoptimize_climbing
        self.NEB_max_images =int(NEB_max_images)

        self.res = Result()
        self.res.new_transition_states = []
    
    def _refineTransitionStates(self, neb, climbing_images):
        """
        refine the transition state candidates.  If at least one is successful
        then return True
        """
        #find the nearest transition state to the transition state candidates
        nrefine = min(self.nrefine_max, len(climbing_images))
        count = 0
        success = False
        for energy, i in climbing_images[:nrefine]:
            count += 1
            print ""
            print "refining transition state from NEB climbing image:", count, "out of", nrefine
            coords = neb.coords[i,:]
            #get guess for initial eigenvector from NEB tangent
            if True:
                eigenvec0 = neb.tangent( neb.energies[i], neb.energies[i-1], neb.energies[i+1],
                                         neb.distance(neb.coords[i,:], neb.coords[i-1,:])[1],
                                         neb.distance(neb.coords[i,:], neb.coords[i+1,:])[1],
                                        )
            
            ret = _refineTS(self.pot, coords, tsSearchParams=self.tsSearchParams, 
                                 eigenvec0=eigenvec0)
            ts_success = ret[0]
            if ts_success:
                #the transition state is good, add it to the list
                tsret, m1ret, m2ret = ret[1:4]
                self.res.new_transition_states.append( (tsret, m1ret, m2ret) )
                success = True
        return success

    def _getNEB(self, *args, **kwargs):
        """
        wrap the actual call to initializing the NEB object so it can be overloaded
        """
        return create_NEB(*args, **kwargs)
#        return NEB(*args, **kwargs)
   
    def _doNEB(self, minNEB1, minNEB2, repetition = 0):
        """
        do NEB between minNEB1 and minNEB2.
        """
        #arrange the coordinates to minimize the distance between them        
        dist, newcoords1, newcoords2 = self.mindist(minNEB1.coords, minNEB2.coords)
        print ""
        
        if repetition == 0: 
            factor = 1.
        else: 
            #change parameters for second repetition
            print "running NEB a second time"
            print "    doubling the number of images"
            print "    doubling the number of steps"
            factor = float(repetition + 1)
        
        print "starting NEB run to try to connect minima", minNEB1._id, minNEB2._id, dist
        neb = self._getNEB(self.pot, newcoords1, newcoords2, 
                         NEBquenchParams=self.NEBquenchParams, 
                         verbose=True, factor=factor, **self.NEBparams)
#        neb.optimize(**NEBquenchParams)
        neb.optimize()
        neb.MakeAllMaximaClimbing()

        if self.reoptimize_climbing > 0:
            print "optimizing climbing images for a small number of steps"
#            NEBquenchParams["nsteps"] = self.reoptimize_climbing
#            neb.optimize(**NEBquenchParams)
            neb.quenchParams["nsteps"] = self.reoptimize_climbing
            neb.optimize()


    
        #get the transition state candidates from the NEB result
        climbing_images = [ (neb.energies[i], i) for i in range(neb.nimages) 
                           if neb.isclimbing[i] ]
        return climbing_images, neb
    
    
    def connect(self, min1, min2):
        """
        1) NEB to find transition state candidates.  
        
        for each transition state candidate:
        
            2) refine the transition state candidates
        
            3) if successful, fall off either side of the transition state
            to find the minima the transition state connects. Add the new 
            transition state and minima to the graph 
        """
        self.NEBattempts = 2
        for repetition in range(self.NEBattempts):
            #do NEB run
            climbing_images, neb = self._doNEB(min1, min2, repetition)
            
            #check results
            nclimbing = len(climbing_images)
            self.res.nclimbing = nclimbing
            print "from NEB search found", nclimbing, "transition state candidates"
            if nclimbing > 0:
                climbing_images = sorted(climbing_images, reverse=True) #highest energies first
                #refine transition state candidates
                self.res.success = self._refineTransitionStates(neb, climbing_images)
            else :
                self.res.success = False
            if self.res.success:
                break
        
        return self.res


#####
#####
#below here only testing routines
#####
#####

def getRandomCoords(pot, natoms):
    import numpy as np
    coords = np.random.uniform(-1,1,natoms*3)*natoms**(1./3)*1.5
    ret = defaults.quenchRoutine(coords, pot.getEnergyGradient)
    return ret

def getPairLJ(natoms=38):
    from pygmin.potentials.lj import LJ
    lj = LJ()
    ret1 = getRandomCoords(lj, natoms)
    ret2 = getRandomCoords(lj, natoms)
    coords1, coords2 = ret1[0], ret2[0]
    E1, E2 = ret1[1], ret2[1]
    
    from pygmin.mindist import MinDistWrapper, minPermDistStochastic
    mindist = MinDistWrapper(minPermDistStochastic, permlist=[range(natoms)])
    mindist(coords1, coords2)
    
    return coords1, coords2, lj, mindist, E1, E2

def test():
    from pygmin.storage import Database
    coords1, coords2, pot, mindist, E1, E2 = getPairLJ()
    db = Database()
    min1 = db.addMinimum(E1, coords1)
    min2 = db.addMinimum(E2, coords2)
    
    
    local_connect = LocalConnect(pot, mindist)
    local_connect.connect(min1, min2)

if __name__ == "__main__":
    test()