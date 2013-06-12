from pele.utils.disconnectivity_graph import DisconnectivityGraph
from pele.storage import Database
from pele.landscape import TSGraph
import pylab as pl
import numpy as np

kbT = 0.75

db = Database(db="tip4p_8.sqlite", createdb=False)

graph = TSGraph(db)

dg = DisconnectivityGraph(graph.graph, db.minima(), subgraph_size=20)
dg.calculate()
dg.plot()

for m in db.minima():
    if m.pgorder != 2:
        print m.pgorder
    m.free_energy = m.energy + kbT * 0.5*m.fvib + kbT*np.log(m.pgorder)

 
for ts in db.transition_states():
#    if ts.pgorder != 2:
    print ts.pgorder
    #assert ts.pgorder == 2    

    ts.free_energy = ts.energy + kbT * 0.5*ts.fvib + kbT*np.log(ts.pgorder) + kbT*np.log(kbT)
    if ts.free_energy > ts.minimum1.free_energy or ts.free_energy > ts.minimum2.free_energy:        
        print "warning, free energy of transition state lower than minimum"
        print ts.free_energy, ts.minimum1.free_energy, ts.minimum2.free_energy
        print ts.energy, ts.minimum1.energy, ts.minimum2.energy
        
    
dg_F = DisconnectivityGraph(graph.graph, db.minima(), energy_attribute="free_energy", subgraph_size=20)
dg_F.calculate()
dg_F.plot()
    
pl.show()

