import numpy as np
from pele.potentials.lj import LJ
from pele.optimize import lbfgs_scipy as quench
from pele.mc import MonteCarlo
from pele.takestep.displace import RandomDisplacement as TakeStep
from pele.takestep.adaptive import AdaptiveStepsize
import sys
from pele.utils.xyz import write_xyz

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
            write_xyz(self.fout, coords, title=str(E))
        self.count += 1

#def center(E, coords, accepted):
 #   c = np.reshape()

natoms = 40

coords = np.random.rand(natoms*3)


lj = LJ() 

ret = quench(coords, lj)
coords = ret.coords

takestep = TakeStep(stepsize=0.1 )

#do equilibration steps, adjusting stepsize
tsadapt = AdaptiveStepsize(takestep)
mc = MonteCarlo(coords, lj, tsadapt, temperature = 0.5)
equilout = open("equilibration","w")
#mc.outstream = equilout
mc.setPrinting(equilout, 1)
mc.run(10000)

#fix stepsize and do production run
mc.takeStep = takestep
mcout = open("mc.out", "w")
mc.setPrinting(mcout, 10)
#print coords from time to time
xyzout = open("out.xyz", "w")
printevent = PrintEvent(xyzout, 300)
mc.addEventAfterStep(printevent.printwrapper)

mc.run(10000)
