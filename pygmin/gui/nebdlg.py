import matplotlib
matplotlib.use("QT4Agg")

from collections import deque
import numpy as np
from PyQt4.QtGui import QDialog, QApplication
import sys

from pygmin.storage import Database
from pygmin.utils.events import Signal
import ui.nebbrowser

def no_event(*args, **kwargs):
    return


# the following lines implement the ability to follow the status of the
# NEB in real time.  the callback is passed to the NEB and plots the
# last several values of the energies of the path.  The implementation
# is quite simple and could easily be a lot faster.
# set follow_neb=False to turn this off
class NEBCallback(object):
    def __init__(self, plw, axes, frq=30, nplots=3):
        self.count = 0
        self.nplots = nplots
        self.data = deque()
        self.frq = frq
        self.process_events = Signal()
        self.plw = plw
        self.axes = axes
        
        
    def __call__(self, energies=None, distances=None, **kwargs):
        self.count += 1
        if self.count % self.frq == 1:
#            print "plotting NEB energies"
            S = np.zeros(energies.shape)
            S[1:] = np.cumsum(distances)
            self.data.append((S, energies.copy()))
            if len(self.data) > self.nplots:
                self.data.popleft()
            self.axes.clear()
            for S, E in self.data:
                line, = self.axes.plot(S, E, "o-")
                # note: if we save the line and use line.set_ydata(E)
                # this would be a lot faster, but we would have to keep
                # track of the y-axis limits manually
#            pl.ylabel("NEB image energy")
#            pl.xlabel("distance along the path")
            self.plw.draw()
            self.process_events()
            import pylab as pl
            pl.pause(.0001)
#            if self.app is not None:
#                self.app.processEvents()



class NEBDialog(QDialog):
    def __init__(self, min1, min2, system):
        super(NEBDialog, self).__init__()
                
        self.ui = ui.nebbrowser.Ui_Form()
        self.ui.setupUi(self)
        self.plw = self.ui.widget
        
        self.system = system
        self.min1 = min1
        self.min2 = min2
        
#        self.minimum_selected = Signal()
        # self.minimum_selected(minim)
    
        self.process_events = Signal()
            

    def do_NEB(self, coords1, coords2):
        throwaway_db = Database()
        min1 = throwaway_db.addMinimum(0., coords1)
        min2 = throwaway_db.addMinimum(1., coords2)
        #use the functions in DoubleEndedConnect to set up the NEB in the proper way
        double_ended = self.system.get_double_ended_connect(min1, min2, 
                                                            throwaway_db, 
                                                            fresh_connect=True)
        local_connect = double_ended._getLocalConnectObject()

        
        
        follow_neb = True
        if follow_neb:
            neb_callback = NEBCallback(self.plw, self.plw.axes)
            neb_callback.process_events.connect(self.process_events)
        else:
            neb_callback = no_event
        
        self.neb =  local_connect._getNEB(self.system.get_potential(),
                                          coords1, coords2, event=neb_callback,
                                          verbose=True,
                                          **local_connect.NEBparams)        
        
        self.neb.optimize()
#        self.nebcoords = self.neb.coords
#        self.nebenergies = self.neb.energies
#        self.ui.oglPath.setCoords(self.neb.coords[0,:], 1)
#        self.ui.oglPath.setCoords(None, 2)
#        self.ui.sliderFrame.setRange(0,self.neb.coords.shape[0]-1)
#        if self.usepymol:
#            self.pymolviewer.update_coords(self.nebcoords, index=1, delete_all=True)

def start():
    print "starting  neb"
    wnd.do_NEB(min1.coords, min2.coords)
    
if __name__ == "__main__":
    from pygmin.systems import LJCluster
    from pygmin.storage import Database
    import pylab as pl
    app = QApplication(sys.argv)
    
    def process_events():
#        print "pausing"
        #pl.pause(.00001)
        app.processEvents()
    
    natoms = 13
    system = LJCluster(natoms)
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = Database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    
    pl.ion()
#    pl.show()
    wnd = NEBDialog(min1, min2, system)    
    wnd.show()
    wnd.process_events.connect(process_events)
    #wnd.process_events = process_events
#    wnd.plw.show()
    
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
        
    sys.exit(app.exec_()) 
        