from pygmin.potentials.maxneib_blj import MaxNeibsBLJ, MaxNeibsBLJSystem
from pygmin.gui import run_gui

natoms = 20
ntypeA = natoms/2
max_neibs = 3
system = MaxNeibsBLJSystem(natoms, ntypeA=ntypeA, max_neibs=max_neibs, rneib=1.7, epsneibs=5.)

dbname = "blj_N%d_NA%d_n%d.db" %(natoms, ntypeA, max_neibs)

run_gui(system, db=dbname)