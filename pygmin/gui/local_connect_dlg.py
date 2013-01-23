import numpy as np
from PyQt4.QtGui import QDialog, QApplication
import sys

from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from pygmin.landscape import Graph
from pygmin.storage import Database
from pygmin.utils.events import Signal
from pygmin.transition_states import findTransitionState, minima_from_ts
from pygmin.landscape.local_connect import _refineTS

import local_connect_browser
from nebdlg import getNEB

class LocalConnectDialog(QDialog):
    def __init__(self, system):
        super(LocalConnectDialog, self).__init__()
        
        self.system = system
#        self.database = database
        
        self.ui = local_connect_browser.Ui_Form()
        self.ui.setupUi(self)
        self.nebwgt = self.ui.wgt_neb
#        self.nebwgt.show()

        self.oglwgt = self.ui.wgt_ogl_slider
        self.oglwgt.setSystem(self.system)
        
        self.nebwgt.on_neb_pick.connect(self.on_neb_pick)


    def _prepareNEB(self, coords1, coords2):
        """setup the NEB object"""
        system = self.system
        throwaway_db = Database()
        min1 = throwaway_db.addMinimum(0., coords1)
        min2 = throwaway_db.addMinimum(1., coords2)
        #use the functions in DoubleEndedConnect to set up the NEB in the proper way
        self.double_ended_connect = system.get_double_ended_connect(min1, min2, 
                                                            throwaway_db, verbosity=0,
                                                            fresh_connect=True)
        self.local_connect = self.double_ended_connect._getLocalConnectObject()
    
        
        
        self.neb =  self.local_connect.create_neb(system.get_potential(),
                                          coords1, coords2,
                                          verbose=True,
                                          **self.local_connect.NEBparams)        
        
#        return neb


    def createNEB(self, coords1, coords2):
        self._prepareNEB(coords1, coords2)
        self.nebwgt.attach_to_NEB(self.neb)
        
        return self.neb


    def on_neb_pick(self, energy, index):
#        print "in local_connect.  you picked E", energy, "index", index
        self.neb_chosen_index = index
        self.oglwgt.setCoordsPath(self.neb.coords, frame=index)

    def runNEB(self):
        self.neb=self.neb.run()
        self.neb_labels = ["NEB path: energy=%f"%(E) for E in self.neb.energies]
        self.show_neb_path()
#        wnd.oglwgt.setCoordsPath(wnd.neb.coords)
    
#    def showFrame(self, i):
#        if hasattr(self, "neb"):
#            self.ui.oglPath.setCoords(self.neb.coords[i,:])
    
#    def on_slider_ogl_sliderMoved(self, index):
#        print "slider moved", index
#        self.oglwgt.setCoords(self.neb.coords[0,:], index=index)

    def on_btn_refineTS_clicked(self):
        self.refine_transition_state()

    def refine_transition_state_old(self):
        print "refining ts"
        coords = self.oglwgt.oglwgt.coords[1].copy()
        tsdata = []
        coordslist = []
        
        def findTS_callback(coords=None, energy=None, rms=None, eigenval=None, **kwargs):
            coordslist.append(coords.copy())
            tsdata.append((energy, rms, eigenval))
        self.ts_result = findTransitionState(coords, self.system.get_potential(),
                                  event=findTS_callback, 
                                  **self.local_connect.tsSearchParams)
        self.ts_coordspath = np.array(coordslist)
        self.ts_labels = ["TS path: energy=%g, rms=%g, eigenval=%g"%(vals) for vals in tsdata]
#        print "tscoords shape", self.ts_coords.shape, len(coordslist)
        self.show_TS_path()
        
        self.pushoff_TS()
    def refine_transition_state(self):
        print "refining ts"
        coords = self.oglwgt.oglwgt.coords[1].copy()
        tsdata = []
        tscoordslist = []
        
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
        pushoff_params["quenchParams"]["events"] = [pushoff_callback]

        
        success, tsret, min1ret, min2ret = _refineTS(self.system.get_potential(), coords,
                    tsSearchParams=tsSearchParams, pushoff_params=pushoff_params)
        
        #process the results of the transition state search
        self.ts_coordspath = np.array(tscoordslist)
        self.ts_labels = ["TS path: energy=%g, rms=%g, eigenval=%g"%(vals) for vals in tsdata]
#        print "tscoords shape", self.ts_coords.shape, len(coordslist)
        self.show_TS_path()
        
        if success:
            #process the results from the pushoff and quenches
            ret1 = min1ret[4]
            ret2 = min2ret[4]

            #the paths from falling off both sides are in pcoordslist.  try to split them up
            i = ret1.nsteps
            self.pushoff_coordspath1 = np.array([tsret.coords.copy()] + pcoordslist[:i])
            self.pushoff_coordspath2 = np.array([tsret.coords.copy()] + pcoordslist[i:])
            data1 = [(tsret.energy, tsret.rms)] + pdata[:i]
            data2 = [(tsret.energy, tsret.rms)] + pdata[i:]
            self.pushoff_labels1 = ["Pushoff left: energy=%g, rms=%g"%(vals) for vals in data1]
            self.pushoff_labels2 = ["Pushoff right: energy=%g, rms=%g"%(vals) for vals in data2]


    def pushoff_TS(self):
        print "falling off either side of the transition state"
        ret = self.ts_result
        #check to make sure it is a valid transition state 
        coords = ret.coords.copy()
        if not ret.success:
            print "transition state search failed"
            return False
            
        if ret.eigenval >= 0.:
            print "warning: transition state has positive lowest eigenvalue", ret.eigenval, ret.energy, ret.rms
            return False
        
        data = []
        coordslist = []
        def pushoff_callback(coords=None, energy=None, rms=None, **kwargs):
            coordslist.append(coords)
            data.append((energy, rms))
        
        kwargs = self.local_connect.pushoff_params.copy()
        kwargs["quenchParams"]["events"] = [pushoff_callback]
        
        #find the minima which this transition state connects
        print "falling off either side of transition state to find new minima"
        ret1, ret2 = minima_from_ts(self.system.get_potential(), coords, n=ret.eigenvec,
            **self.local_connect.pushoff_params)
        #get the Results objects
        ret1 = ret1[4] 
        ret2 = ret2[4]


        #the paths from falling off both sides are in coordslist.  try to split them up
        i = ret1.nsteps
        self.pushoff_coordspath1 = np.array([ret.coords.copy()] + coordslist[:i])
        self.pushoff_coordspath2 = np.array([ret.coords.copy()] + coordslist[i:])
        data1 = [(ret.energy, ret.rms)] + data[:i]
        data2 = [(ret.energy, ret.rms)] + data[i:]
        self.pushoff_labels1 = ["Pushoff left: energy=%g, rms=%g"%(vals) for vals in data1]
        self.pushoff_labels2 = ["Pushoff right: energy=%g, rms=%g"%(vals) for vals in data2]



#        self.pushoff_coordspath = np.array(coordslist)
#        self.pushoff_labels = ["Pushoff path: energy=%g, rms=%g"%(vals) for vals in data]

    def show_pushoff_path1(self):
        self.oglwgt.setCoordsPath(self.pushoff_coordspath1, labels=self.pushoff_labels1)
    def show_pushoff_path2(self):
        self.oglwgt.setCoordsPath(self.pushoff_coordspath2, labels=self.pushoff_labels2)
    
    def show_neb_path(self, frame=0):
        self.oglwgt.setCoordsPath(self.neb.coords, frame=frame, labels=self.neb_labels)

    def show_TS_path(self):
        self.oglwgt.setCoordsPath(self.ts_coordspath, labels=self.ts_labels)
        
        
        

def start():
    print "starting  neb"
    wnd.createNEB(x1, x2)
    wnd.runNEB()
    
    
if __name__ == "__main__":
    from pygmin.systems import LJCluster
    from nebdlg import getNEB
    from OpenGL.GLUT import glutInit

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
    wnd = LocalConnectDialog(system)   
    wnd.show()
    wnd.nebwgt.process_events.connect(process_events)
    

    glutInit()

    #initilize the NEB and run it.
    #we have to do it through QTimer because the gui has to 
    #be intitialized first... I don't really understand it 
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)

    sys.exit(app.exec_()) 
