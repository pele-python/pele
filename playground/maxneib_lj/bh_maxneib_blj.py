from past.utils import old_div
from pele.potentials.maxneib_blj import MaxNeibsBLJ, MaxNeibsBLJSystem
from pele.gui import run_gui
import numpy as np


natoms = 50
ntypeA = int(old_div(natoms, 2))
max_neibs = 7
only_AB_neibs = True
rneib = 0.74 * (3.0 / 4 * np.pi)

periodic = True
if periodic:
    rho = 1.0
    boxl = (float(natoms) / rho) ** (1.0 / 3)
    print(boxl)
else:
    boxl = None


system = MaxNeibsBLJSystem(
    natoms,
    ntypeA=ntypeA,
    max_neibs=max_neibs,
    neib_crossover=0.3,
    rneib=rneib,
    epsneibs=6.0,
    epsAB=1.0,
    epsB=0.01,
    epsA=0.01,
    sigA=1.3,
    sigB=1.3,
    sigAB=1.0,
    boxl=boxl,
    only_AB_neibs=only_AB_neibs,
)
# system.params.basinhopping.outstream = None

if only_AB_neibs:
    onlyAB = "_onlyAB"
else:
    onlyAB = ""
if periodic:
    textboxl = "_boxl%.2f" % boxl
else:
    textboxl = ""
dbname = "blj_N%d_NA%d_n%d%s%s_rneib%.2f_newsig.db" % (
    natoms,
    ntypeA,
    max_neibs,
    textboxl,
    onlyAB,
    rneib,
)
print(dbname)

gui = True
if gui:
    run_gui(system, db=dbname)
else:
    db = system.create_database(dbname)

    bh = system.get_basinhopping(database=db)

    bh.run(1000000)
