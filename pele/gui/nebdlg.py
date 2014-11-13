import matplotlib
matplotlib.use("QT4Agg")

from collections import deque
import numpy as np
from PyQt4.QtGui import QDialog, QApplication, QWidget, QVBoxLayout
import sys
from itertools import izip
from matplotlib.lines import Line2D
from matplotlib.patches import Circle


from pele.storage import Database
from pele.utils.events import Signal
import ui.nebbrowser

def no_event(*args, **kwargs):
    return


# the following lines implement the ability to follow the status of the
# NEB in real time.  the callback is passed to the NEB and plots the
# last several values of the energies of the path.  The implementation
# is quite simple and could easily be a lot faster.
# set follow_neb=False to turn this off
class NEBCallback(object):
    def __init__(self, plw, axes, neb=None, frq=30, nplots=3):
        self.count = 0
        self.nplots = nplots
        self.data = deque()
        self.frq = frq
        self.process_events = Signal()
        self.plw = plw
        self.axes = axes
        self.neb = neb
        
        self.on_coords_select = Signal()

        self.points_list = []
        def on_pick_tmp(event):
            self.on_pick(event)
            return
        self.plw.mpl_connect('pick_event', on_pick_tmp)
     
    def on_pick(self, event):
        for points, nebdata in izip(self.points_list, self.data):
            if event.artist == points:
                #nebdata is in the form  (S, E, stepnum)
                S, E, stepnum = nebdata
                ind = event.ind[0]
                print "you picked a point with energy", E[ind], "index", ind
                self.on_coords_select(energy=E[ind], index=ind) 
        
    def __call__(self, energies=None, distances=None, stepnum=None, **kwargs):
        self.count += 1
        if self.count % self.frq == 1:
#            print "plotting NEB energies"
            S = np.zeros(energies.shape)
            S[1:] = np.cumsum(distances)
            self.data.append((S, energies.copy(), stepnum))
            if len(self.data) > self.nplots:
                self.data.popleft()
            self.axes.clear()
            
            
            self.points_list = []
            for S, E, N in self.data:
                line, = self.axes.plot(S, E, "-", label=str(N))
                points = self.axes.scatter(S, E, picker=5)
                self.points_list.append(points) 
                # note: if we save the line and use line.set_ydata(E)
                # this would be a lot faster, but we would have to keep
                # track of the y-axis limits manually
            self.points_list = [self.points_list[-1]] #only enable interactive for the last line plotted

            self.axes.set_ylabel("NEB image energy")
            self.axes.set_xlabel("distance along the path")
            self.axes.legend(title="step num")
            leg = self.axes.legend(title="step num", fancybox='True')
            leg.get_frame().set_alpha(0.5)
            self.plw.draw()
            self.process_events()




class NEBWidget(QWidget):
    def __init__(self, parent=None):
        super(NEBWidget, self).__init__(parent=parent)
                
        self.ui = ui.nebbrowser.Ui_Form()
        self.ui.setupUi(self)
        self.plw = self.ui.widget
        
#        self.plw.axes.set_ylabel("NEB image energy")
#        pl.xlabel("distance along the path")
        
#        self.system = system
#        self.min1 = min1
#        self.min2 = min2
        
#        self.minimum_selected = Signal()
        # self.minimum_selected(minim)
    
        #the function which tells the eventloop to continue processing events
        #if this is not set, the plots will stall and not show until the end of the NEB run.
        #ideally it is the function app.processEvents where app is returned by
        #app = QApplication(sys.argv)
        self.process_events = Signal()
        
        self.on_neb_pick = Signal()
        self.on_neb_pick.connect(self.on_pick)
    
    def highlight2(self, x, y):
        """draw a circle around x, y"""
#        if hasattr(self, "highlight_circle"):
        try:
            self.highlight_circle.remove()
        except Exception: pass
        self.highlight_circle=Circle((x, y),radius=0.5, fill=False)
        self.plw.axes.add_patch(self.highlight_circle)
        
    
    def highlight1(self, x, y):
        ylim = self.plw.axes.get_ylim()
#        if hasattr(self, "highlight_line"):
        try:
            self.highlight_line.remove()
        except Exception: pass
        self.highlight_line = Line2D([x,x],list(ylim), ls='--', c='k')
#        self.highlight_line.set_data([x,x],list(ylim))
                                                        
        self.plw.axes.add_line(self.highlight_line)
#        self.plw.axes.plot([x,x], list(ylim), 'k')
    
    def highlight(self, index, style=1):
        """draw a vertical line to highlight a particular point in the neb curve"""
        if index < 0:
            #I would like to delete the line here, but I don't know how to do it easily
            return
        from matplotlib.lines import Line2D
        S, energies, tmp = self.neb_callback.data[-1]
        x = S[index]
        y = energies[index]
        if style == 2:
            self.highlight2(x, y)
        else:
            self.highlight1(x, y)
        
        self.plw.draw()
        self.process_events()
        
    
    def on_pick(self, index=None, *args, **kwargs):
        self.highlight(index)

    def attach_to_NEB(self, neb):
        neb_callback = NEBCallback(self.plw, self.plw.axes)
        self.neb_callback = neb_callback
        neb_callback.process_events.connect(self.process_events)
        neb_callback.on_coords_select.connect(self.on_neb_pick)      
        neb.update_event.connect(neb_callback)

class NEBDialog(QDialog):
    def __init__(self, *args, **kwargs):
        super(NEBDialog, self).__init__(*args, **kwargs)
        self.nebwgt = NEBWidget(parent=self)

        vbox = QVBoxLayout()
        vbox.addWidget(self.nebwgt)
#        vbox.addWidget(self.mpl_toolbar)
        self.setLayout(vbox)



def getNEB(coords1, coords2, system):
    """setup the NEB object"""
    throwaway_db = Database()
    min1 = throwaway_db.addMinimum(0., coords1)
    min2 = throwaway_db.addMinimum(1., coords2)
    #use the functions in DoubleEndedConnect to set up the NEB in the proper way
    double_ended = system.get_double_ended_connect(min1, min2, 
                                                        throwaway_db, 
                                                        fresh_connect=True)
    local_connect = double_ended._getLocalConnectObject()

    
    
    neb =  local_connect._getNEB(system.get_potential(),
                                      coords1, coords2,
                                      verbose=True,
                                      **local_connect.NEBparams)        
    
    return neb


def start():
    print "starting  neb"
    neb = getNEB(x1, x2, system)
#    wnd.do_NEB(min1.coords, min2.coords)
    wnd.attach_to_NEB(neb)
    neb.optimize()
    
if __name__ == "__main__":
    from pele.systems import LJCluster
    from pele.storage import Database
    import pylab as pl
    app = QApplication(sys.argv)
    
    def process_events():
        app.processEvents()
    
    #setup system
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = Database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    #setup neb dialog
    pl.ion()
#    pl.show()
    dlg = NEBDialog()
    wnd = dlg.nebwgt   
    dlg.show()
    wnd.process_events.connect(process_events)

    #initilize the NEB and run it.
    #we have to do it through QTimer because the gui has to 
    #be intitialized first... I don't really understand it 
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
        
    sys.exit(app.exec_()) 
        