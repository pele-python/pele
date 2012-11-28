import numpy as np

from pygmin.transition_states import InterpolatedPathDensity

__all__ = ["smoothPath"]

def smoothPath(path, mindist, density=5.):
    fullpath = []
    coords1 = path[0] 
    for i in range(1,len(path)):
        #coords1 = path[i-1]
        coords2 = path[i]
        dist, coords1, coords2 = mindist(coords1, coords2)
        newpath = InterpolatedPathDensity(coords1, coords2, dist, density)
        fullpath += list(newpath)
        #print "dist", dist, len(newpath), len(fullpath)
        
        coords1 = coords2
    return fullpath
    