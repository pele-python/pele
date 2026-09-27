from past.utils import old_div
from pele.potentials.maxneib_blj import MaxNeibsBLJ, MaxNeibsBLJSystem
from pele.gui import run_gui

natoms = 20
ntypeA = old_div(natoms, 2)
max_neibs = 4.5
system = MaxNeibsBLJSystem(
    natoms,
    ntypeA=ntypeA,
    max_neibs=max_neibs,
    rneib=1.7,
    epsneibs=12.0,
    epsAB=1.0,
    epsB=0.01,
    epsA=0.01,
    neib_crossover=0.05,
)

dbname = "blj_N%d_NA%d_n%d.db" % (natoms, ntypeA, max_neibs)

run_gui(system)  # , db=dbname)
