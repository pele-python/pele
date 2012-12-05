import numpy as np
from pygmin.potentials import LJ
from pygmin.utils import xyz
from pygmin.angleaxis import rigidbody

# read in coordinates from xyz file
ref = xyz.read_xyz(open("water.xyz"))
xyz.write_xyz(open("test.xyz", "w"), coords = ref.coords)
# lookup table for atom masses
mass_lookup = {'O': 16., 'H': 1.}

#ref.coords[:] *= 3

# now define a new rigid body system
rb_sites = []
for atomtype, x, i in zip(ref.atomtypes, ref.coords, xrange(len(ref.atomtypes))):
    # every 3rd atom, define a new reigid molecule
    if i%3 == 0:
        rb = rigidbody.RigidFragment()
        rb_sites.append(rb)                
    rb.add_atom(atomtype, x, mass_lookup[atomtype])

# finalize the rigid body setup
for rb in rb_sites:
    rb.finalize_setup()

# define a new rigid body system    
rbsystem = rigidbody.RBSystem()
rbsystem.add_sites(rb_sites)

print len(rbsystem.sites), len(rbsystem.indices)

rbcoords = rbsystem.coords_adapter(np.zeros(len(rbsystem.sites)*6))

for site, com in zip(rbsystem.sites, rbcoords.posRigid):
    com[:] = ref.coords[site.indices[0]] - site.atom_positions[0]
     
# for simplicity just use a lj potential here
pot = LJ()

# get the flattened coordinate array
print pot.getEnergy(ref.coords.flatten())
rbpot = rigidbody.RBPotentialWrapper(rbsystem, pot)
print rbpot.getEnergy(rbcoords.coords)
#print rbpot.getEnergyGradient(rbcoords.coords)

