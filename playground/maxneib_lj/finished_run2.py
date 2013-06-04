from pele.potentials.maxneib_blj import MaxNeibsBLJ, MaxNeibsBLJSystem
from pele.gui import run_gui
import numpy as np


natoms = 50
ntypeA = int(natoms/2)
max_neibs = 3
only_AB_neibs = True
rneib = 1.3 #.74 * (3./4*np.pi)

periodic = False
if periodic:
    rho = 1.
    boxl = (float(natoms) / rho)**(1./3)
    print boxl
else:
    boxl=None



system = MaxNeibsBLJSystem(natoms, ntypeA=ntypeA, max_neibs=max_neibs, 
                           neib_crossover=.3,
                           rneib=rneib,
                           epsneibs=6., epsAB=1., epsB=0.01, epsA=0.01, 
                           sigB=1.,
#                           sigA=1.3, sigB=1.3, sigAB=1.,
                           boxl=boxl,
                           only_AB_neibs=only_AB_neibs)
#system.params.basinhopping.outstream = None

if only_AB_neibs:
    onlyAB = "_onlyAB"
else:
    onlyAB = ""
if periodic:
    textboxl = "_boxl%.2f" % boxl
else:
    textboxl = ""
dbname = "blj_N%d_NA%d_n%d%s%s_rneib%.2f.db" %(natoms, ntypeA, max_neibs, textboxl, onlyAB, rneib)
print dbname

gui = True
if gui:
    run_gui(system, db=dbname)
else:
    db=system.create_database(dbname)
    
    bh = system.get_basinhopping(database=db)
    
    bh.run(1000000)
