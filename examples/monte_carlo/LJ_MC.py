import numpy as np
from potentials.lj import LJ
from optimize.quench import quench
from mc import MonteCarlo
from take_step.random_displacement import takeStep
import sys
from printing.print_atoms_xyz import printAtomsXYZ as printxyz

class PrintEvent(object):
    def __init__(self, fout, printfrq = 1):
        self.fout = fout
        self.printfrq = printfrq
        self.count = 0
        
    def center(self, coords):
        natoms = len(coords)/3
        for i in range(3):
            com = np.sum(coords[i::3])/natoms
            coords[i::3] -= com

        
    def printwrapper(self, E, coords, accepted):
        if self.count % self.printfrq == 0: 
            self.center(coords)
            printxyz(self.fout, coords, line2=str(E))
        self.count += 1

#def center(E, coords, accepted):
 #   c = np.reshape()

natoms = 40

coords = np.random.rand(natoms*3)


lj = LJ() 

ret = quench(coords, lj.getEnergyGradient)
coords = ret[0]

takestep = takeStep(stepsize=0.1 )

#do equilibration steps, adjusting stepsize
takestep.useAdaptiveStep()
mc = MonteCarlo(coords, lj, takestep, temperature = 0.5)
equilout = open("equilibration","w")
#mc.outstream = equilout
mc.setPrinting(equilout, 1)
mc.run(10000)

#fix stepsize and do production run
takestep.useFixedStep()
mcout = open("mc.out", "w")
mc.setPrinting(mcout, 10)
#print coords from time to time
xyzout = open("out.xyz", "w")
printevent = PrintEvent(xyzout, 300)
mc.addEventAfterStep(printevent.printwrapper)

mc.run(10000)