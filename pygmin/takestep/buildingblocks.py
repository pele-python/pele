import numpy as np
from pygmin.utils import rotations
from pygmin.utils import vec3

def uniform_displace(stepsize, coords, indices=None):
    if(indices):
        for i in indices:
            coords[i] += stepsize*rotations.vec_random()
        return
    
    for x in coords:
        x += stepsize*rotations.vec_random()
        
def rotate(stepsize, coords, indices=None):
    if(indices):
        for i in indices:
            rotations.takestep_aa( coords[i], stepsize )
        return
    
    for x in coords:
        rotations.takestep_aa( x, stepsize )
        
def reduced_coordinates_displace(stepsize, lattice_matrix, coords, indices=None):
    ilattice = vec3.invert3x3(lattice_matrix) # inverse_lattice
    if(indices):
        for i in indices:
            coords[i] += np.dot(ilattice, stepsize*rotations.vec_random())
        return
            
    for x in coords:
        x += np.dot(ilattice, stepsize*rotations.vec_random())
   
    