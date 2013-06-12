import numpy as np
import sandbox_wrapper as sandbox

nmol = 3
nsites = 3*3
natoms = nmol*2 # num. com coords + num. aa coords
r = 5
s = .5
os= .5
myas = .5

sandbox.setup_commons(nmol, nsites, r, s, os, myas)

sandbox.input("data")

com = np.array([ 0.,0.,0., 1.,0.,0., 0.,1.,0. ])
aacoords = np.array([ 1.,0.,0., 0.,1.,0., 0.,0.,1.,] )
coords = np.zeros( 3*2*nmol, np.float64)
coords[0:3*nmol] = com[:]
coords[3*nmol:] = aacoords[:]

printlist = []
printlist.append( coords.copy() )
grad, E = sandbox.getenergygradient(coords, [2*nmol])
print E, grad


sandbox.takestep(coords, False, True, s, os, myas,[natoms])
printlist.append( coords.copy() )

sandbox.takestep(coords, True, False, s, os, myas,[natoms])
printlist.append( coords.copy() )




from printing.print_atoms_xyz import printAtomsXYZ as printxyz
with open("out.xyz", "w") as fout:
    for x in printlist:
        printxyz(fout, x)
        
        
        
        
        