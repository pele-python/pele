import Pyro4

nruns = 5
print "I will do %d connect runs before exiting"%nruns

print "reading uri to connect to pyros server"
fl=open("pyros.uri")
uri = fl.read()
fl.close()

print "get connect manager"
connect_manager=Pyro4.Proxy(uri) 

system = connect_manager.get_system()
pot = system.get_potential()

# forwarding new minima to server
def minimum_added(minimum):
    minid = connect_manager.add_minimum(minimum.energy, minimum.coords)
    # store the id of minimum on serverside for quick access later on
    minimum.parentid = minid

# forwarding new ts to server 
def ts_added(ts):
    id1 = ts.minimum1.parentid
    id2 = ts.minimum2.parentid
    
    tsid = connect_manager.add_ts(id1, id2, ts.energy, ts.coords)

# create a fake database in memory
db = system.create_database(db=":memory:")

# connect to events and forward them to server
db.on_minimum_added.connect(minimum_added)
db.on_ts_added.connect(ts_added)
    
for i in xrange(nruns):
    print "Obtain a new job"
    id1, coords1, id2, coords2 = connect_manager.get_connect_job()
    
    print "add minima to local database"
    min1 = db.addMinimum(pot.getEnergy(coords1), coords1)
    min2 = db.addMinimum(pot.getEnergy(coords2), coords2)
    
    print "run double ended connect"
    connect = system.get_double_ended_connect(min1, min2, db, fresh_connect=True)
    connect.connect()

print "finished sucessfully!"
print "Data collected during run:"
print db.session.qery(Minimum).count(), "minima"
print db.session.qery(Minimum).count(), "transition states"
