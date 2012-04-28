import numpy as np
import sandbox_wrapper as sandbox

nmol = 3
nsites = 3*3
natoms = nmol*2 # num. com coords + num. aa coords
r = 5
s = .5
os= .5
myas = .5

sandbox.sandbox_wrapper.setup_commons(nmol, nsites, r, s, os, myas)

sandbox.sandbox_wrapper.input("data")

com = np.array([ 0.,0.,0., 1.,0.,0., 0.,1.,0. ])
aacoords = np.array([ 1.,0.,0., 0.,1.,0., 0.,0.,1.,] )
coords = np.zeros( 3*2*nmol, np.float64)
coords[0:3*nmol] = com[:]
coords[3*nmol:] = aacoords[:]
grad, E = sandbox.sandbox_wrapper.getenergygradient(coords, [2*nmol])
print E, grad
