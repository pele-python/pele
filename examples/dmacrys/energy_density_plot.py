'''
Created on Jul 9, 2012

@author: vr274
'''

from pygmin.utils.rbtools import *
from pygmin.utils import lattice
import pylab as pl
import numpy as np

def energy_density_plot(storage):
    rho=[]
    ener=[]
    for minimum in storage.data:        
        ca = CoordsAdapter(nrigid=2, nlattice=6, coords=minimum.coords)
        ml = lattice.lowerTriangular(ca.lattice)
        volume = np.dot(ml[:,0],np.cross(ml[:,1], ml[:,2]))
        rho.append(1./volume)
        ener.append(minimum.energy)        
    return rho, ener           
        
if __name__ == '__main__':
    import pickle
    save = pickle.load(open("storage"))
    rho, ener = energy_density_plot(save)
    print ener    
    pl.plot(rho, ener, 'x')
    pl.axis(ymax = -0.7)
    pl.show() 
    
