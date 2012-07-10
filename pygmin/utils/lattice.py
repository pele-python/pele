import numpy as np
def lowerTriangular(coords):
    return np.array(        
        [[coords[0],    0.,        0.],
        [ coords[1], coords[3],    0.],
        [ coords[2], coords[4], coords[5]]])

def volume(coords):
    ml = lowerTriangular(coords)
    return np.dot(ml[:,0],np.cross(ml[:,1], ml[:,2]))

    