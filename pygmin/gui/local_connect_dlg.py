import numpy as np
from PyQt4.QtGui import QDialog, QApplication
import sys

from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from pygmin.landscape import Graph
from pygmin.storage import Database
from pygmin.utils.events import Signal

import local_connect_browser
from nebdlg import getNEB

class LocalConnectDialog(QDialog):
    def __init__(self, system, database):
        super(LocalConnectDialog, self).__init__()
        
        self.system = system
        self.database = database
        
        self.ui = local_connect_browser.Ui_Form()
        self.ui.setupUi(self)
        self.nebwgt = self.ui.wgt_neb
#        self.nebwgt.show()

        self.oglwgt = self.ui.wgt_ogl_slider
        self.oglwgt.setSystem(self.system)

    def createNEB(self, coords1, coords2):
        self.neb = getNEB(coords1, coords2, self.system)
        self.nebwgt.attach_to_NEB(self.neb)
        
        return self.neb


    def on_neb_pick(self, energy, index):
#        print "in local_connect.  you picked E", energy, "index", index
        self.oglwgt.setCoordsPath(self.neb.coords, frame=index)

    def runNEB(self):
        self.neb.optimize()
        wnd.oglwgt.setCoordsPath(wnd.neb.coords)

    
    def showFrame(self, i):
        if hasattr(self, "neb"):
            self.ui.oglPath.setCoords(self.neb.coords[i,:])
    
#    def on_slider_ogl_sliderMoved(self, index):
#        print "slider moved", index
#        self.oglwgt.setCoords(self.neb.coords[0,:], index=index)


def start():
    print "starting  neb"
    wnd.createNEB(x1, x2)
#    neb = getNEB(x1, x2, system)
##    wnd.do_NEB(min1.coords, min2.coords)
#    wnd.nebwgt.attach_to_NEB(neb)

#    print "running neb"
#    neb.optimize()
    wnd.runNEB()
    
    
if __name__ == "__main__":
    from pygmin.systems import LJCluster
    from nebdlg import getNEB
    import pylab as pl
    from OpenGL.GLUT import glutInit

    app = QApplication(sys.argv)
    
    def process_events():
        app.processEvents()
    
    #setup system
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 3.
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = Database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    #setup neb dialog
    pl.ion()
#    pl.show()
    wnd = LocalConnectDialog(system, db)   
    wnd.show()
#    wnd.nebwgt.show()
    wnd.nebwgt.process_events.connect(process_events)
    
    def on_pick(*args, **kwargs):
        return wnd.on_neb_pick(*args, **kwargs)
    wnd.nebwgt.on_neb_pick.connect(on_pick)

    glutInit()

    #initilize the NEB and run it.
    #we have to do it through QTimer because the gui has to 
    #be intitialized first... I don't really understand it 
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)

    sys.exit(app.exec_()) 
