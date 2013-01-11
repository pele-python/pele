"""
BH and connections using OpenMM Amber Potential       

26-dec-2012 -- all working here 
 
todo: integrate with GUI in guisystem.py  
"""
import openMMamberPot as amb  
import pygmin.utils.amber as readAmb

from amberSystem import AMBERsystem  

from pygmin.systems import BaseSystem
from pygmin.mindist import minPermDistStochastic, MinDistWrapper, ExactMatchCluster
from pygmin.transition_states import orthogopt

from simtk.unit import  angstrom 
    
'''  ------------------------------------------------------------------- '''
import numpy as np
def run_basinhopping(x0):
    
    natoms = len(coords)/3 
#    sysAmb = MySystem(natoms)
    sysAmb = AMBERsystem(natoms)
        
#    database = sysAmb.create_database()
    
    import pygmin.gui.run as gr
    gr.run_gui(sysAmb)

    
    from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
    takeStepRnd   = RandomDisplacement( stepsize=2 )
    tsAdaptive = AdaptiveStepsizeTemperature(takeStepRnd, interval=10, verbose=True)

    sysAmb.params.database.accuracy = 1e-3
    sysAmb.params.basinhopping["temperature"] = 10.0        
#    sysAmb.params.basinhopping["insert_rejected"] = True      
    # todo - how do you save N lowest?    

#    bh = sysAmb.get_basinhopping(database=database, takestep = tsAdaptive, coords=x0)    
    bh = sysAmb.get_basinhopping(database=database, takestep = takeStepRnd, coords=x0)    
          
    # todo -- how to set run parameters like tolerances; 
    #         have hardcoded insert_rejected in basinhopping.py 
    
    print 'Running BH .. '
    bh.run(5)
    print "Number of minima found = ", len(database.minima())
    min0 = database.minima()[0]
    print "lowest minimum found has energy = ", min0.energy
    return sysAmb, database



'''  ------------------------------------------------------------------- '''

def run_double_ended_connect(sys, database):
    #connect the all minima to the lowest minimum
    minima = database.minima()
    min1 = minima[0]
    
#    for min2 in minima[1:]:
    for min2 in minima[1:]:
        connect = sys.get_double_ended_connect(min1, min2, database)
        connect.connect()

'''  ------------------------------------------------------------------- '''
from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from pygmin.landscape import Graph
import matplotlib.pyplot as plt

def make_disconnectivity_graph(database):
    graph = Graph(database).graph
    dg = DisconnectivityGraph(graph, nlevels=3, center_gmin=True)
    dg.calculate()
    dg.plot()
    plt.show()

'''
  -------------------------------------------------------------------
'''
def test_potential(coords):
    pot = amb.amberPot()       
    e = pot.getEnergy(coords)
    print 'Energy (kJ/mol) = '
    print e
    
    e, g = pot.getEnergyGradient(coords)
    gnum = pot.NumericalDerivative(coords, eps=1e-6)

    print 'Energy (kJ/mol) = '
    print e 
    print 'Analytic Gradient = '
    print g[60:65] 
    print 'Numerical Gradient = '
    print gnum[60:65] 

    print 'Num vs Analytic Gradient =' 
    print np.max(np.abs(gnum-g)), np.max(np.abs(gnum))
    print np.max(np.abs(gnum-g)) / np.max(np.abs(gnum))

'''
  ==============================================================
'''
if __name__ == "__main__":
    
    from pygmin.gui import run as gr    

    # create new amber system 
    sysAmb   = AMBERsystem('coords.prmtop', 'coords.inpcrd')        
    # create new database  
    database = sysAmb.create_database()
    # run gui 
    gr.run_gui(sysAmb)
    
    exit()

    # read one conformation from pdb file 
    from simtk.openmm.app import pdbfile as openmmpdb
    pdb = openmmpdb.PDBFile('coords.pdb')
    
    openmmpdb.PDBFile.writeFile(   topology, positions, file, modelIndex)
    
    coords = pdb.getPositions() / angstrom   
    coords = np.reshape(np.transpose(coords), 3*len(coords), 1)

    test_potential(coords)
        
    # Basin Hopping     
    mysys, database = run_basinhopping(coords)
    
    # Connect run
    # run_double_ended_connect(mysys, database)
    
    # Disconn 
    # make_disconnectivity_graph(database)

#    import pygmin.gui.run as gr
#    gr.run_gui(MySystem)
    
    
    
    
