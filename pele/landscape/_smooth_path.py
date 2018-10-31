from pele.transition_states import InterpolatedPathDensity

__all__ = ["smoothPath", "smooth_path"]


def smooth_path(path, mindist, density=5., interpolator=None):
    """return a smooth (interpolated) path
    
    useful for making movies.  Especially for min-ts-min pathways
    returned by DoubleEndedConnect.
    
    Parameters
    ----------
    path : list of arrays
        the pathway to smooth
    mindist : callable
        function that returns the mindist distance between two structures ::
        
            dist, newcoords1, newcoords2 = mindist(coords1, coords2)
    
    density : float, optional
        how dense to do the smoothing.
    interpolator : callable, optional
        allows to specify a custom interpolation routine (e.g. for angle axis)
        
    """
    fullpath = []
    coords1 = path[0]
    for i in range(1, len(path)):
        coords2 = path[i]
        dist, coords1, coords2 = mindist(coords1, coords2)
        newpath = InterpolatedPathDensity(coords1, coords2, dist, density, interpolator=interpolator)
        # skip the first image because it is already in fullpath 
        fullpath += list(newpath)[1:]

        coords1 = coords2
    return fullpath

smoothPath = smooth_path # for backwards compatibility
    
