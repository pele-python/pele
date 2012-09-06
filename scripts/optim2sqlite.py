from pygmin.storage.database import Database, Minimum, TransitionState
import time
db = Database(db="storage.sqlite")
 
i=1
minima={}
t0 = time.time()
tt = t0

print "Reading minima"
for line in open("min.data"):
    energy, frequency, pgorder, itx, ity, itz = line.split()
    min1 = Minimum(float(energy), None)
    db.session.add(min1)
    minima[i]=min1
    i+=1
print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
tt = time.time()


print "Reading transition states"
for line in open("ts.data"):
    energy, frequency, pgorder, min1, min2, itx, ity, itz = line.split()
    ts = TransitionState(float(energy), None, minima[int(min1)], minima[int(min2)])
    db.session.add(ts)
    #db.addTransitionState(float(energy), None, minima[int(min1)], minima[int(min2)], commit=False)
print "%.1f seconds"%(time.time() - tt)
tt = time.time()
print "Commiting changes to database"
db.session.commit()
print "%.1f seconds"%(time.time() - tt)
print "Done after %.1f seconds"%(time.time() - t0)

from  pygmin.NEB.graph import Graph
print "generating graph"
tt = time.time()
graph = Graph(db)
print "%.1f seconds"%(time.time() - tt)

#print "printing minima"
#for i in db.minima():
#    print i