import numpy as np
from PyQt4.QtGui import QDialog, QApplication, QListWidgetItem
from PyQt4 import QtCore
import sys

from pele.storage import Database
from pele.landscape.local_connect import _refineTS

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class _TransitionStateView(object):
    """this hold all the data necessary for ConnectExplorer to save transition state search data"""
    def __init__(self, nebindex, ts_coordspath, ts_labels, pushoff_coordspath, pushoff_labels):
        self.nebindex = nebindex
        self.ts_coordspath = ts_coordspath
        self.ts_labels = ts_labels
        self.pushoff_coordspath = pushoff_coordspath
        self.pushoff_labels = pushoff_labels

class _TSListItem(QListWidgetItem):
    def __init__(self, nebindex, *args, **kwargs):
        text="nebindex %d"% nebindex
        QListWidgetItem.__init__(self, text)
        
        self.tsview = _TransitionStateView(nebindex, *args, **kwargs)
        
        


class ConnectExplorerDialog(QDialog):
    def __init__(self, system, app, parent=None):
        super(ConnectExplorerDialog, self).__init__(parent=parent)
        
        self.system = system
        self.app = app
#        self.database = database
        import connect_explorer_ui

        self.ui = connect_explorer_ui.Ui_Form()
        self.ui.setupUi(self)
        self.nebwgt = self.ui.wgt_neb
#        self.nebwgt.show()

        self.oglwgt = self.ui.wgt_ogl_slider
        self.oglwgt.setSystem(self.system)
        
        self.ts_list = self.ui.list_ts
        
        self.nebwgt.on_neb_pick.connect(self.on_neb_pick)
        
        QtCore.QObject.connect(self.oglwgt.slider, QtCore.SIGNAL(_fromUtf8("sliderMoved(int)")),
                               self.highlight_frame)

        self.oglview = "None"
    
    def reset(self):
        """clear everything and start again"""
        self.ts_list.clear()
    
    def highlight_frame(self, index):
        if self.oglview == "neb":
            self.nebwgt.highlight_frame(index)
        else:
            self.nebwgt.highlight_frame(-1)

    def on_neb_pick(self, index):
#        print "in local_connect.  you picked E", energy, "index", index
        self.neb_chosen_index = index
        self.show_neb_path(frame=index)
        self.ui.wgt_neb.highlight_frame(index)

    def set_nebrunner(self, nebrunner):
        self.ui.wgt_neb.update_gui(nebrunner)
        self.nebrunner = nebrunner
        self.neb_labels = ["energy = %f"%e for e in nebrunner.energies[-1]]
        self.show_neb_path()
        
    def on_refine_all_ts(self):
        print "refining all"
        self.nebrunner.neb.neb.MakeAllMaximaClimbing()
        climbing_images = [i for i in range(len(self.nebrunner.path)) if self.nebrunner.neb.neb.isclimbing[i]]
#        print "climbing images", climbing
        for nebindex in climbing_images:
            self.refine_transition_state(nebindex=nebindex)            

    def on_refine_transition_state(self):
        self.refine_transition_state()

    def refine_transition_state(self, nebindex=None):
        print "refining ts"
        if nebindex is None:
            #figure out which image to start from
            if self.oglview != "neb":
                print "choose which NEB image to start from"
                return
    #            raise Exception("choose which NEB image to start from")
            nebindex = self.oglwgt.get_slider_index()
            
        self.highlight_frame(-1)
        coords = self.nebrunner.path[nebindex].copy()
        self.ui.wgt_neb.highlight_point(nebindex)
        self.app.processEvents()
    
        tsdata = []
        tscoordslist = []
        
        self.local_connect = self.nebrunner.local_connect
        #setup the callback function for findTransitionState
        def findTS_callback(coords=None, energy=None, rms=None, eigenval=None, **kwargs):
            tscoordslist.append(coords.copy())
            tsdata.append((energy, rms, eigenval))
            
        tsSearchParams = self.local_connect.tsSearchParams.copy()
        tsSearchParams["event"] = findTS_callback
        
        #setup the callback function for the pushoff
        pdata = []
        pcoordslist = []
        def pushoff_callback(coords=None, energy=None, rms=None, **kwargs):
            pcoordslist.append(coords)
            pdata.append((energy, rms))
        pushoff_params = self.local_connect.pushoff_params.copy()
        quencher = self.system.get_minimizer(events=[pushoff_callback])
        pushoff_params["quench"] = quencher

        
        success, tsret, min1ret, min2ret = _refineTS(self.system.get_potential(), coords,
                    tsSearchParams=tsSearchParams, pushoff_params=pushoff_params)
        
        #process the results of the transition state search
        self.ts_coordspath = np.array(tscoordslist)
        self.ts_labels = ["TS path: energy=%g, rms=%g, eigenval=%g"% vals for vals in tsdata]
#        print "tscoords shape", self.ts_coords.shape, len(coordslist)
        self.show_TS_path()
        
        if success:
            #process the results from the pushoff and quenches
            ret1 = min1ret[4]
            ret2 = min2ret[4]

            #the paths from falling off both sides are in pcoordslist.  try to split them up
            i = ret1.nsteps
            pushoff_coordspath1 = [tsret.coords.copy()] + pcoordslist[:i]
            pushoff_coordspath2 = [tsret.coords.copy()] + pcoordslist[i:]
            data1 = [(tsret.energy, tsret.rms)] + pdata[:i]
            data2 = [(tsret.energy, tsret.rms)] + pdata[i:]
            pushoff_labels1 = ["Pushoff left: energy=%g, rms=%g"% vals for vals in data1]
            pushoff_labels2 = ["Pushoff right: energy=%g, rms=%g"% vals for vals in data2]
            #combine them together with one in reversed order.
            self.pushoff_coordspath = list(reversed(pushoff_coordspath1)) + pushoff_coordspath2
            self.pushoff_labels = list(reversed(pushoff_labels1)) + pushoff_labels2
            self.pushoff_coordspath = np.array(self.pushoff_coordspath)
            
            #make a
            tsitem = _TSListItem(nebindex, self.ts_coordspath, self.ts_labels, self.pushoff_coordspath, self.pushoff_labels)
            self.ts_list.addItem(tsitem)
        

    def load_ts_view(self, tsitem):
        #highlight the neb image we started from
        ts = tsitem.tsview
        self.ui.wgt_neb.highlight_point(ts.nebindex)
        
        self.ts_coordspath = ts.ts_coordspath
        self.ts_labels = ts.ts_labels
        self.pushoff_coordspath = ts.pushoff_coordspath
        self.pushoff_labels = ts.pushoff_labels
        
        self.show_TS_path()
        

            
    def on_list_ts_selected(self, item):
        self.load_ts_view(item)   



    def show_pushoff_path(self):
        self.oglwgt.setCoordsPath(self.pushoff_coordspath, labels=self.pushoff_labels)
        self.oglview = "pushoff"
#    def show_pushoff_path2(self):
#        self.oglwgt.setCoordsPath(self.pushoff_coordspath2, labels=self.pushoff_labels2)
#        self.oglview = "pushoff1"
    
    def show_neb_path(self, frame=0):
        self.oglwgt.setCoordsPath(self.nebrunner.path, frame=frame, labels=self.neb_labels)
        self.oglview = "neb"
        

    def show_TS_path(self):
        self.oglwgt.setCoordsPath(self.ts_coordspath, labels=self.ts_labels, frame=-1)
        self.oglview = "ts"
        
        
        

def start():
    print "starting  neb"
    from neb_explorer import NEBRunner
    runner = NEBRunner(app, system)
    runner.run(x1, x2)
    wnd.set_nebrunner(runner)
    
    
if __name__ == "__main__":
    from pele.systems import LJCluster
    from nebdlg import getNEB
    from OpenGL.GLUT import glutInit

    app = QApplication(sys.argv)
    
    def process_events():
        app.processEvents()
    
    #setup system
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 10.
    system.params.double_ended_connect.local_connect_params.NEBparams.image_density = 3.
#    system.params.double_ended_connect.local_connect_params.NEBparams.adaptive_nimages = 5.
    system.params.double_ended_connect.local_connect_params.NEBparams.reinterpolate = 400
    system.params.double_ended_connect.local_connect_params.NEBparams.max_images = 40
    x1, e1 = system.get_random_minimized_configuration()[:2]
    x2, e2 = system.get_random_minimized_configuration()[:2]
    db = Database()
    min1 = db.addMinimum(e1, x1)
    min2 = db.addMinimum(e2, x2)
    
    #setup neb dialog
    wnd = ConnectExplorerDialog(system, app)   
    wnd.show()

    glutInit()

    #initilize the NEB and run it.
    #we have to do it through QTimer because the gui has to 
    #be intitialized first... I don't really understand it 
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)

    sys.exit(app.exec_()) 
