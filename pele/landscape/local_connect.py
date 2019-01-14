import logging

from pele.optimize import Result
from pele.transition_states import findTransitionState, minima_from_ts
from pele.transition_states import NEBDriver

__all__ = ["LocalConnect"]

logger = logging.getLogger("pele.connect")


def _refineTS(pot, coords, tsSearchParams=None, eigenvec0=None, pushoff_params=None):
    """
    find nearest transition state to NEB climbing image.  Then fall
    off the transition state to find the associated minima.
    
    This would naturally be a part of DoubleEndedConnect.  I separated it
    to make it more easily parallelizable.      
    """
    if pushoff_params is None: pushoff_params = dict()
    if tsSearchParams is None: tsSearchParams = dict()
    # run ts search algorithm
    kwargs = dict(tsSearchParams.items())
    ret = findTransitionState(coords, pot, eigenvec0=eigenvec0, **kwargs)

    # check to make sure it is a valid transition state 
    coords = ret.coords
    if not ret.success:
        logger.info("transition state search failed")
        return False, ret, None, None

    if ret.eigenval >= 0.:
        logger.info("transition state has positive lowest eigenvalue, skipping: %s %s %s", ret.eigenval, ret.energy,
                    ret.rms)
        logger.info("         not adding transition state")
        return False, ret, None, None

    # find the minima which this transition state connects
    logger.info("falling off either side of transition state to find new minima")
    ret1, ret2 = minima_from_ts(pot, coords, n=ret.eigenvec,
                                **pushoff_params)

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
    nrefine_max : int
        the maximum number of NEB transition state candidates to refine
    reoptimize_climbing : int
        the number of iterations to use for re-optimizing the climbing images
        after the NEB is done.
    pushoff_params : int
        parameters for detemining how to find the minima on either side of 
        a transition state
    verbosity : int
        this controls how many status messages are printed.  (not really
        implemented yet)
    
    
    Notes
    -----
    this class takes two minima as input, does an NEB run to find transition state
    candidates, then refines those candidates into transition states.  Finally, 
    we fall off either side of the transition states to fine the minima on either
    side.
    
    This is the core routine of DoubleEndedConnect.  It is separated out in order
    to make parallelization easier.  This class intentionally has no knowledge of the
    global landscape (database, graph, etc.).
    
    See Also
    --------
    DoubleEndedConnect : the routine from which local connect is generally called
    pele.transition_states.NEB : one of the core routines
    pele.transition_states.NEBDriver : the wrapper which sets up NEB
    pele.transition_states.findTransitionState : one of the core routine
    
    """

    def __init__(self, pot, mindist, tsSearchParams=None, verbosity=1, NEBparams=None, nrefine_max=100,
                 reoptimize_climbing=0, pushoff_params=None, create_neb=NEBDriver):
        if pushoff_params is None: pushoff_params = dict()
        if NEBparams is None: NEBparams = dict()
        if tsSearchParams is None: tsSearchParams = dict()
        self.pot = pot
        self.mindist = mindist
        self.tsSearchParams = tsSearchParams
        self.verbosity = int(verbosity)
        self.nrefine_max = nrefine_max

        self.NEBparams = NEBparams
        self.reoptimize_climbing = reoptimize_climbing

        self.pushoff_params = pushoff_params

        self.res = Result()
        self.res.new_transition_states = []
        self.create_neb = create_neb

    def _refineTransitionStates(self, neb, climbing_images):
        """
        refine the transition state candidates.  If at least one is successful
        then return True
        """
        # find the nearest transition state to the transition state candidates
        nrefine = min(self.nrefine_max, len(climbing_images))
        count = 0
        success = False
        for energy, i in climbing_images[:nrefine]:
            count += 1
            logger.info("")
            logger.info("refining transition state from NEB climbing image: %s %s %s", count, "out of", nrefine)
            coords = neb.coords[i, :]
            # get guess for initial eigenvector from NEB tangent
            if True:
                eigenvec0 = neb.tangent(neb.energies[i], neb.energies[i - 1], neb.energies[i + 1],
                                        neb.distance(neb.coords[i, :], neb.coords[i - 1, :])[1],
                                        neb.distance(neb.coords[i, :], neb.coords[i + 1, :])[1],
                )

            ret = _refineTS(self.pot, coords, tsSearchParams=self.tsSearchParams,
                            eigenvec0=eigenvec0, pushoff_params=self.pushoff_params)
            ts_success = ret[0]
            if ts_success:
                # the transition state is good, add it to the list
                tsret, m1ret, m2ret = ret[1:4]
                self.res.new_transition_states.append((tsret, m1ret, m2ret))
                success = True
        return success

    def _doNEB(self, minNEB1, minNEB2, repetition=0):
        """
        do NEB between minNEB1 and minNEB2.
        """
        # arrange the coordinates to minimize the distance between them        
        dist, newcoords1, newcoords2 = self.mindist(minNEB1.coords, minNEB2.coords)
        logger.info("")

        if repetition == 0:
            factor = 1.
        else:
            # change parameters for second repetition
            logger.info("running NEB a second time")
            logger.info("    doubling the number of images")
            logger.info("    doubling the number of steps")
            factor = float(repetition + 1)

        logger.info("starting NEB run to try to connect minima %s %s %s", minNEB1.id(), minNEB2.id(), dist)

        neb = self.create_neb(self.pot, newcoords1, newcoords2,
                              factor=factor, **self.NEBparams)
        neb = neb.run()

        neb.MakeAllMaximaClimbing()

        if self.reoptimize_climbing > 0:
            logger.info("optimizing climbing images for a small number of steps")
            # NEBquenchParams["nsteps"] = self.reoptimize_climbing
            # neb.optimize(**NEBquenchParams)
            neb.quenchParams["nsteps"] = self.reoptimize_climbing
            neb.optimize()



        # get the transition state candidates from the NEB result
        climbing_images = [(neb.energies[i], i) for i in range(neb.nimages)
                           if neb.isclimbing[i]]
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
            # do NEB run
            climbing_images, neb = self._doNEB(min1, min2, repetition)
            self.neb = neb

            # check results
            nclimbing = len(climbing_images)
            self.res.nclimbing = nclimbing
            logger.info("from NEB search found %s %s", nclimbing, "transition state candidates")

            # if the NEB path has maxima, refine those into transition states
            if nclimbing > 0:
                climbing_images = sorted(climbing_images, reverse=True)  # highest energies first
                # refine transition state candidates
                self.res.success = self._refineTransitionStates(neb, climbing_images)
            else:
                self.res.success = False
            if self.res.success:
                break

        return self.res


# ####
# ####
# below here only testing routines
#####
#####

def getPairLJ(natoms=38):  # pragma: no cover
    from pele.systems import LJCluster

    system = LJCluster(natoms)
    ret1 = system.get_random_minimized_configuration()
    ret2 = system.get_random_minimized_configuration()
    coords1, coords2 = ret1[0], ret2[0]
    E1, E2 = ret1[1], ret2[1]

    mindist = system.get_mindist()
    mindist(coords1, coords2)

    return coords1, coords2, system.get_potential(), mindist, E1, E2


def test():  # pragma: no cover
    from pele.storage import Database

    coords1, coords2, pot, mindist, E1, E2 = getPairLJ()
    db = Database()
    min1 = db.addMinimum(E1, coords1)
    min2 = db.addMinimum(E2, coords2)

    local_connect = LocalConnect(pot, mindist)
    local_connect.connect(min1, min2)


if __name__ == "__main__":
    #    logger.basicConfig(level=logger.DEBUG)
    test()

