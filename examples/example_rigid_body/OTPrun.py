import numpy as np
import potentials.rigid_bodies.molecule as molecule
import potentials.rigid_bodies.sandbox as sandbox
from potentials.rigid_bodies.take_step import RBTakeStep
import rotations as rot
import copy


nmol = 5

#define the molecule types.
#here use only one type, LWOTP
otp = molecule.setupLWOTP()


# define the interaction matrix for the system.
# for LWOTP there is only one atom type, so this is trivial
from potentials.lj import LJ
lj = LJ()
interaction_matrix = [[lj]]


#set up a list of molecules
mols = [otp for i in range(nmol)]

#set up the RBSandbox object
mysys = sandbox.RBSandbox(mols, interaction_matrix)
nsites = mysys.nsites

#get an initial set of coordinates
coords = np.zeros(2*3*nmol, np.float64)
coords[0:nmol*3] = np.random.uniform(-1,1,[nmol*3]) * 1.3*(nsites)**(1./3)
for i in range(nmol):
    k = nmol*3 + 3*i
    coords[k : k + 3] = rot.random_aa()


#set up the takestep routine
takestep = RBTakeStep()

#set up the class to save lowest energy structures
from storage.savenlowest import SaveN
saveit = SaveN(100)

#set up basinhopping
from basinhopping import BasinHopping
bh = BasinHopping(coords, mysys, takestep.take_step, storage=saveit.insert )

#run basin hopping
bh.run(40)



#print the saved coords
fname = "otp.xyz"
print "saving xyz coords to", fname
from printing.print_atoms_xyz import printAtomsXYZ as printxyz
with open(fname, "w") as fout:
    for E, coords in saveit.data:
        xyz = mysys.getxyz(coords)
        printxyz( fout, xyz, atom_type=["N", "O", "O"])


