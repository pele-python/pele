'''
Created on Aug 1, 2012

@author: vr274
'''

from pygmin.storage.database import Storage

if __name__ == '__main__':
    db = Storage(db="../pygmin/NEB/test.db")
    
    minindex={}
    out = open("min.data", "w")
    i=1
    print "foo"
    for m in db.minima():
        minindex[m]=i
        i+=1
        out.write("%f 0.0 1 0.0 0.0 0.0\n"%(m.energy))
        print "%f 0.0 1 0.0 0.0 0.0"%(m.energy)
        
    out = open("ts.data", "w")
    for ts in db.transition_states():
        out.write("%f 0.0 1 %d %d 0.0 0.0 0.0\n"%(ts.energy, minindex[ts.minimum1], minindex[ts.minimum2]))
        print "%f 0.0 1 %d %d 0.0 0.0 0.0"%(ts.energy, minindex[ts.minimum1], minindex[ts.minimum2])
        