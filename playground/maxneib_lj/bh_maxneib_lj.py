from pele.potentials.maxneib_lj import MaxNeibsLJ, MaxNeibsLJSystem

natoms = 20
max_neibs = 3
system = MaxNeibsLJSystem(natoms, max_neibs=max_neibs, rneib=1.7, epsneibs=5.)
#system.params.basinhopping.outstream = None

dbname = "lj_N%d_n%d.db" %(natoms, max_neibs)

db=system.create_database(dbname)

bh = system.get_basinhopping(database=db)

bh.run(1000000)
