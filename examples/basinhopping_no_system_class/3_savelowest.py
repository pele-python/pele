# ###########################################################
# Example 3: Saving the coordinates as an xyz file
# ###########################################################
import numpy as np
import pele.potentials.lj as lj
import pele.basinhopping as bh
from pele.takestep import displace
from pele.storage.savenlowest import SaveN as saveit


natoms = 12
coords = np.random.random(3 * natoms)

potential = lj.LJ()
step = displace.RandomDisplacement(stepsize=0.5)

storage = saveit(nsave=10)
opt = bh.BasinHopping(coords, potential, step, storage=storage.insert)
opt.run(100)

with open("lowest", "w") as fout:
    fout.write(str(natoms) + "\n")
    fout.write(str(storage.data[0].energy) + "\n")
    atoms = storage.data[0].coords.reshape(natoms, 3)
    for a in atoms:
        fout.write("LA " + str(a[0]) + " " + str(a[1]) + " " + str(a[2]) + "\n")

############################################################
# some visualization
############################################################
try:
    import pele.utils.pymolwrapper as pym

    pym.start()
    frame = 1
    print storage.data
    for minimum in storage.data:
        coords = minimum.coords.reshape(natoms, 3)
        pym.draw_spheres(coords, "A", frame)
        frame += 1
except:
    print "Could not draw using pymol, skipping this step" 