import numpy as np
from pele.angleaxis._otp_cluster import OTPCluster


from pele.angleaxis.tests.test_otp import _x1, _x2
from pele.transition_states import NEBDriver

def getdb(nmol=20):
    system = OTPCluster(nmol)
    db = system.create_database("test.sqlte")
    bh = system.get_basinhopping(db)
    bh.run(20)
    

def do_mindist(nmol=20):
    system = OTPCluster(nmol)
    db = system.create_database("test.sqlte")
    mindist = system.get_mindist()
    m0 = db.minima()[0]
    for m1 in db.minima()[1:6]:
        d, x1, x2 = mindist(m0.coords, m1.coords)
        print d

def do_NEB(nmol=20):
    system = OTPCluster(nmol)
    db = system.create_database("test.sqlte")
    mindist = system.get_mindist()
    pot = system.get_potential()
    m1 = db.minima()[0]
    for m2 in db.minima()[1:3]:
        neb = NEBDriver(pot, m1.coords, m2.coords)
        neb.run()

def do_connect(nmol=20):
    system = OTPCluster(nmol)
    db = system.create_database("test.sqlte")
    m1 = db.minima()[0]
    for m2 in db.minima()[1:5]:
        connect = system.get_double_ended_connect(m1, m2, db, fresh_connect=True)
        connect.connect()

def do_ts_search(nmol=20):
    from pele.transition_states import findTransitionState
    system = OTPCluster(nmol)
    db = system.create_database("test.sqlte")
    orthogopt = system.get_orthogonalize_to_zero_eigenvectors()
    for ts in db.transition_states():
        coords = db.transition_states()[0].coords.copy() 
        coords += np.random.uniform(-.5,.5, coords.size)
    
        print system.params.double_ended_connect.local_connect_params.tsSearchParams
        findTransitionState(coords, system.get_potential(), orthogZeroEigs=orthogopt, 
                            **system.params.double_ended_connect.local_connect_params.tsSearchParams)
    

if __name__ == "__main__":
#     getdb()
#     do_mindist()    
    do_connect()
#     do_ts_search()