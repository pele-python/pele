import numpy as np

from pygmin.transition_states import NEB, NEBPar, InterpolatedPath

__all__ = ["create_NEB"]

def create_NEB(pot, coords1, coords2, image_density=10, max_images=40,
                iter_density=15, 
                NEBquenchParams=dict(), 
                verbose=False, factor=1, parallel=False, ncores=4, **NEBparams):
    """
    a wrapper function to do the interpolation and set up the nudged elastic band object
    
    Parameters
    ----------
    image_density : float
        how many NEB images per unit distance to use.
    coords1, coords2 : array
        the structures to connect with the band
    max_images : int
        the maximum number of NEB images
    iter_density : float
    NEBquenchParams : dict
        parameters passed to the NEB minimization routine
    NEBparams : dict
        NEB setup parameters.  These are passed directly to the NEB class.  Use
        NEBquenchParams for parameters related to the optimization of the band.
    verbose : bool
    factor : int
        The number of images is multiplied by this factor.  If the number of 
        images is already at it's maximum, then the number of iterations is 
        multiplied by this factor instead
    parallel : bool
        if True, then use class NEBPar to evaluate the image potentials in parallel
    ncores : int
        the number of cores to use.  Ignored if parallel is False
    
    Returns
    -------
    neb : an initialized NEB object
    
    See Also
    ---------
    NEB
    InterpolatedPath
    """
    #determine the number of images to use
    dist = np.linalg.norm(coords1 - coords2)
    nimages = int(max(1., dist) * image_density * factor)
    if nimages > max_images:
        nimages = max_images

    #determine the number of iterations
    NEBquenchParams = NEBquenchParams.copy()
    if NEBquenchParams.has_key("nsteps"):
        niter = NEBquenchParams["nsteps"]
    else:
        niter = int(iter_density * nimages)
        NEBquenchParams["nsteps"] = niter

    #if nimages is already max_images then increasing the number
    #of images with factor will have no effect.  so double the number of steps instead
    if factor > 1. and nimages == max_images:
        niter *= factor
        NEBquenchParams["nsteps"] = niter

    
    if verbose:    
        print "    NEB: nimages", nimages
        print "    NEB: nsteps ", niter

    if parallel:
        return NEBPar(InterpolatedPath(coords1, coords2, nimages), 
                   pot, quenchParams=NEBquenchParams, ncores=ncores, **NEBparams)
    else:
        return NEB(InterpolatedPath(coords1, coords2, nimages), 
                   pot, quenchParams=NEBquenchParams, **NEBparams)

