from pygmin.potentials.maxneib_blj import MaxNeibsBLJ, MaxNeibsBLJSystem

natoms = 20
ntypeA = natoms/2
max_neibs = 3
system = MaxNeibsBLJSystem(natoms, ntypeA=ntypeA, max_neibs=max_neibs, rneib=1.7, epsneibs=5.)
#system.params.basinhopping.outstream = None

dbname = "blj_N%d_NA%d_n%d.db" %(natoms, ntypeA, max_neibs)

db=system.create_database(dbname)

bh = system.get_basinhopping(database=db)

bh.run(1000000)
