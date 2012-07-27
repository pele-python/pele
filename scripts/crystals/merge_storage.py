'''
Created on Jul 27, 2012

@author: vr274
'''
from pygmin.utils.rbtools import *
from pygmin.storage import savenlowest
from pygmin.utils import crystals
import pickle
    
def compareMinima(min1, min2):
    from pygmin.utils import rbtools
    ca1 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min1.coords)
    ca2 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min2.coords)
    match = crystals.compareStructures(ca1, ca2)
    if not match:
        print "Found minimum with similar energy but different structure"
    return match
    

save = savenlowest.SaveN(nsave=1000, accuracy=1e-3)
save.compareMinima = compareMinima

import sys
for i in sys.argv[1:]:
    print i
    save2 = pickle.load(open(i, "r"))
    for m in save2.data:
        save.insert(m.E, m.coords)

pickle.dump(save, open("storage", "w"))
