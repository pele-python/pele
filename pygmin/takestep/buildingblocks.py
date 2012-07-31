import numpy as np
from pygmin.utils import rotations
from pygmin.utils import vec3

def uniform_displace(stepsize, coords):
    for x in coords:
        x += stepsize*rotations.vec_random()
        
def rotate(stepsize, coords):
    for x in coords:
        rotations.takestep_aa( x, stepsize )
        
def reduced_coordinates_displace(stepsize, lattice_matrix, coords):
    ilattice = vec3.invert3x3(lattice_matrix) # inverse_lattice
    for x in coords:
        x += np.dot(ilattice, stepsize*rotations.vec_random())
   
    