import numpy as np

from pele.storage import Database
from playground.xy_model.xy_model_system import XYModlelSystem
from pele.gui.run import run_gui

dbname="xy_10x10.sqlite"

def create_system(L, dbname):
    print "testing whether", dbname, "exists"
    try:
        # if the database already exists get the phases
        db = Database(dbname, createdb=False)
        print dbname, "exists.  getting phases"
        phases = db.get_property("phases").value()
        dim = db.get_property("dim").value()
        assert dim[0] == L
    except IOError:
        print dbname, "doesn't exist, generating random phases"
        phases=None
        dim=[L,L]
    system = XYModlelSystem(dim=dim, phi_disorder=np.pi, phases=phases)
    return system

def run_gui_disorder(L=10, dbname=None):
    if dbname is None:
        dbname = "xy_%dx%d.sqlite" %(L,L)
    system = create_system(L=L, dbname=dbname)
    db = system.create_database(dbname)
    run_gui(system, db=db)

run_gui_disorder(L=12)