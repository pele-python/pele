import matplotlib
matplotlib.use("QT4Agg")
import traceback    
import sys
import copy
import numpy as np

from PyQt4 import QtCore, QtGui, Qt

from pygmin.gui.MainWindow import Ui_MainWindow 
from pygmin.gui.bhrunner import BHRunner
from pygmin.landscape import Graph
from pygmin.gui.dlg_params import DlgParams
from pygmin.config import config
from pygmin.gui.ui.dgraph_dlg import DGraphDialog
#from pygmin.gui.connect_explorer_dlg import ConnectExplorerDialog
from pygmin.gui.connect_run_dlg import ConnectViewer
from pygmin.gui.takestep_explorer import TakestepExplorer
from pygmin.gui.normalmode_browser import NormalmodeBrowser


def excepthook(ex_type, ex_value, traceback_obj):
    """ redirected exception handler """
    
    errorbox = QtGui.QMessageBox()
    msg = "An unhandled exception occurred:\n"+str(ex_type) + "\n\n"\
                     + str(ex_value) + "\n\nTraceback:\n----------"
    for line in traceback.format_tb(traceback_obj):
        msg += "\n" + line
    errorbox.setText(msg)
    errorbox.setStandardButtons(QtGui.QMessageBox.Ignore | QtGui.QMessageBox.Cancel)
    errorbox.setDefaultButton(QtGui.QMessageBox.Cancel)
    if errorbox.exec_() == QtGui.QMessageBox.Cancel:
        raise

class MySelection(object):
    """keep track of which minima have been selected and whether those coordinates have been modified
    """
    def __init__(self):
        self.minimum1 = None
        self.minimum1 = None
        self.coords1 = None
        self.coords1 = None

class QTSInList(QtGui.QListWidgetItem):
    def __init__(self, ts):
        text="%.4f (%d<-%d->%d)"%(ts.energy, ts._minimum1_id, ts._id, ts._minimum2_id)
        QtGui.QListWidgetItem.__init__(self, text)
    
        self.coords = ts.coords
        self.ts = ts
        self.tsid = id(ts)
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.ts.energy > item2.ts.energy

class MinimumStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of minima in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, minimum):
        text="%.4f (%d)"%(minimum.energy, minimum._id)
        super(MinimumStandardItem, self).__init__(text)
        self.minimum = minimum
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.minimum.energy < item2.minimum.energy

class TransitionStateStandardItem(Qt.QStandardItem):
    """defines an item to populate the lists of minima in the GUI
    
    These items will be collected in a QStandardItemModel and viewed using
    QListView
    """
    def __init__(self, ts):
        text="%.4f (%d<-%d->%d)"%(ts.energy, ts._minimum1_id, ts._id, ts._minimum2_id)
        super(TransitionStateStandardItem, self).__init__(text)
        self.ts = ts
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.ts.energy < item2.ts.energy

    def __getattr__(self, name):
        """return the transition state if the minimum is asked for"""
        if name == "minimum":
            return self.ts
        return super(TransitionStateStandardItem, self).__getattr__(name)

class MinimumStandardItemModel(Qt.QStandardItemModel):
    """a class to manage the list of minima for display in the gui"""
    def __init__(self, nmax=None, **kwargs):
        super(MinimumStandardItemModel, self).__init__(**kwargs)
        self.nmax = nmax # the maximum number of minima
        self.issued_warning = False

    def set_nmax(self, nmax):
        self.nmax = nmax
          
    def appendRow(self, item, *args, **kwargs):
        Qt.QStandardItemModel.appendRow(self, item, *args, **kwargs)
        if self.nmax is not None:
            nrows = self.rowCount()
            if nrows > self.nmax:
                if not self.issued_warning:
                    print "warning: limiting the number of minima displayed in the gui to", self.nmax
                    self.issued_warning = True
                # choose an item to remove from the list.  we can't usume it's totally sorted
                # because it might have been a while since it was last sorted
                candidates = [(self.item(r).minimum.energy, r) for r in xrange(nrows-10,nrows)]
                toremove = max(candidates)
                self.takeRow(toremove[1])
                
                

class MainGUI(QtGui.QMainWindow):
    """
    this is the main class for the pygmin gui
    
    Parameters
    ----------
    app : 
        the application object returned by QtGui.QApplication()
    systemtype : system class object
        the system class
    """
    def __init__(self, app, systemtype, parent=None):
        QtGui.QWidget.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.systemtype = systemtype
        self.transition=None
        self.app = app
        self.double_ended_connect_runs = []
        self.pick_count = 0
        

        # set up the sorting of minima
        self.need_sorting = False
        self.need_sorting_minima = False
        self.need_sorting_ts = False
        self._sort_timer = QtCore.QTimer()
        self._sort_timer.timeout.connect(self._delayed_sort)
        
        self.minima_selection = MySelection()
        self.ts_selected = None
        
        # set up the lists of minima in the gui 
        self.minima_list_model = MinimumStandardItemModel()
        self.ui.list_minima_main.setModel(self.minima_list_model)
#        self.ui.list_minima_main.show()
        self.ui.listMinima1.setModel(self.minima_list_model)
        self.ui.listMinima2.setModel(self.minima_list_model)
        
        # set up the list of transition states in the gui
        self.ts_list_model = MinimumStandardItemModel()
        self.ui.list_TS.setModel(self.ts_list_model)
        
        self.NewSystem()

        # determine the maximum number of minima to keep in the lists.
        # this must be done after NewSystem() is called
        try:
            nmax = self.system.params.gui.list_nmax
            self.minima_list_model.set_nmax(nmax)
            self.ts_list_model.set_nmax(nmax)
        except AttributeError:
            pass

        #try to load the pymol viewer.
        try:
            self.usepymol = self.system.params.gui.use_pymol
        except (KeyError, AttributeError):
            self.usepymol = config.getboolean("gui", "use_pymol")
        if self.usepymol:
            try:
                from pymol_viewer import PymolViewer
                self.pymolviewer = PymolViewer(self.system.load_coords_pymol)
            except (ImportError or NotImplementedError):
                self.usepymol = False
            
        if self.usepymol == False:
            #note: glutInit() must be called exactly once.  pymol calls it
            #during pymol.finish_launching(), so if we call it again it will
            #give an error. On the other hand, if we're not using pymol we 
            #must call it.
            from OpenGL.GLUT import glutInit
            glutInit()
        
    def NewSystem(self):
        """
        this is called to initialize the system with a database
        """
        self.system = self.systemtype()
        self.connect_db()
#        db = self.system.create_database()
#        self.system.database = db
#        self.system.database.on_minimum_added.connect(self.NewMinimum)
#        self.system.database.on_minimum_removed.connect(self.RemoveMinimum)
#        self.system.database.on_ts_removed.connect(self.RemoveTS)
#        self.system.database.on_ts_added.connect(self.NewTS)
#        self.minima_list_model.clear()
#        self.ui.list_TS.clear()
            
    def on_action_edit_params_triggered(self, checked=None):
        if checked is None: return
        self.paramsdlg = DlgParams(self.system.params)
        self.paramsdlg.show()
        
    def processEvents(self):
        self.app.processEvents()
     
    def on_action_db_connect_triggered(self, checked=None):
        """
        launch a file browser to connect to an existing database
        """
        if checked is None: return
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Open File', '.')
        self.connect_db(filename)

    def connect_db(self, filename=":memory:"):
        """
        connect to an existing database at location filename
        """
        self.minima_list_model.clear()
        self.ts_list_model.clear()

        db = self.system.create_database(db=filename)
        self.system.database = db
        #add minima to listWidged.  do sorting after all minima are added
        for minimum in self.system.database.minima():
            self.NewMinimum(minimum, sort_items=False)
        self._sort_minima()
        self.NewTS(self.system.database.transition_states())

        self.system.database.on_minimum_added.connect(self.NewMinimum)
        self.system.database.on_minimum_removed(self.RemoveMinimum)
        self.system.database.on_ts_added.connect(self.NewTS)
        self.system.database.on_ts_removed.connect(self.RemoveTS)
    
    def on_list_minima_main_clicked(self, index):
        item = self.minima_list_model.item(index.row())
        minimum = item.minimum
        self.SelectMinimum(minimum)

    def SelectMinimum(self, minimum):
        """when you click on a minimum in the basinhopping tab
        """
        print "selecting minimum", minimum._id, minimum.energy
        self.ui.ogl_main.setSystem(self.system)
        self.ui.ogl_main.setCoords(minimum.coords)
        self.ui.ogl_main.setMinimum(minimum)
        self.ui.oglTS.setSystem(self.system)
        #self.ui.oglTS.setCoords(item.coords)
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=1, delete_all=True)

    def on_listMinima1_clicked(self, index):
        item = self.minima_list_model.item(index.row())
        minimum = item.minimum
        self._SelectMinimum1(minimum)

    def on_listMinima2_clicked(self, index):
        item = self.minima_list_model.item(index.row())
        minimum = item.minimum
        self._SelectMinimum2(minimum)

    def _SelectMinimum1(self, minimum):
        """by minimum"""
        print "selecting minimum 1", minimum._id, minimum.energy
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=1)
#        self.ui.oglPath.setMinimum(minimum, index=1)
        self.minima_selection.minimum1 = minimum
        self.minima_selection.coords1 = minimum.coords
        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=1)

    def _SelectMinimum2(self, minimum):
        """by minimum"""
        print "selecting minimum 2", minimum._id, minimum.energy
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=2)
#        self.ui.oglPath.setMinimum(minimum, index=2)
        self.minima_selection.minimum2 = minimum
        self.minima_selection.coords2 = minimum.coords

        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=2)
    
    def get_selected_minima(self):
        """return the two minima that have been chosen in the gui"""
        m1, m2 = self.minima_selection.minimum1, self.minima_selection.minimum2
        if m1 is None or m2 is None:
            raise Exception("you must select two minima first")
        return m1, m2
 
    def get_selected_coords(self):
        """return the two sets of coordinates that have been chosen in the gui
        
        note, that these may not be the same as what is stored in the minimum.  E.g. they
        may be the aligned structures
        """
        coords1, coords2 = self.minima_selection.coords1, self.minima_selection.coords2
        if coords1 is None or coords2 is None:
            raise Exception("you must select two minima first")
        return coords1, coords2

#    def on_list_TS_currentItemChanged(self, item):
#        self.show_TS(item.ts)
    def on_list_TS_clicked(self, index):
        item = self.ts_list_model.item(index.row())
        ts = item.ts
        self.ts_selected = ts
        self.show_TS(ts)

    def show_TS(self, ts):
        """
        show the transition state and the associated minima in the 3d viewer
        """
        self.ui.oglTS.setSystem(self.system)
        m1 = ts.minimum1
        m2 = ts.minimum2
        #  put them in best alignment
        mindist = self.system.get_mindist()
        dist, m1coords, tscoords = mindist(m1.coords, ts.coords)
        dist, m2coords, tscoords = mindist(m2.coords, ts.coords)
        self.tscoordspath = np.array([m1coords, tscoords, m2coords])
        labels  = ["minimum: energy " + str(m1.energy) + " id " + str(m1._id)]
        labels += ["ts: energy " + str(ts.energy)]
        labels += ["minimum: energy " + str(m2.energy) + " id " + str(m2._id)]
        self.ui.oglTS.setCoordsPath(self.tscoordspath, frame=1, labels=labels)
        if self.usepymol:
            self.pymolviewer.update_coords(self.tscoordspath, index=1, delete_all=True)


    def on_btnAlign_clicked(self, clicked=None):
        """use mindist to align the minima.
        
        called when the align button is pressed
        """
        if clicked is None: return
        coords1, coords2 = self.get_selected_coords()
        align = self.system.get_mindist()
        pot = self.system.get_potential()
        print "energy before alignment", pot.getEnergy(coords1), pot.getEnergy(coords2)
        dist, coords1, coords2 = align(coords1, coords2)
        print "energy after alignment", pot.getEnergy(coords1), pot.getEnergy(coords2)
        print "best alignment distance", dist
        self.ui.oglPath.setCoords(coords1, index=1)
        self.ui.oglPath.setCoords(coords2, index=2)
        self.minima_selection.coords1 = coords1
        self.minima_selection.coords2 = coords2
        if self.usepymol:
            self.pymolviewer.update_coords([coords1], index=1)
            self.pymolviewer.update_coords([coords2], index=2)
    
    def on_btnNEB_clicked(self, clicked=None):
        """do an NEB run (not a connect run).  Don't find best alignment first"""
        if clicked is None: return
        coords1, coords2 = self.get_selected_coords()
        from neb_explorer import NEBExplorer
        
        if not hasattr(self, "nebexplorer"):
            self.nebexplorer = NEBExplorer(system=self.system, app=self.app, parent=self)
        self.nebexplorer.show()
        self.nebexplorer.new_neb(coords1, coords2)
        
#    def showFrame(self, i):
#        if hasattr(self, "nebcoords"):
#            self.ui.oglPath.setCoords(self.nebcoords[i,:])
    
    def on_minimum_picked(self, min1):
        """called when a minimimum is clicked on in the graph or disconnectivity graph"""
        if (self.pick_count % 2) == 0:
            self._SelectMinimum1(min1)
        else:
            self._SelectMinimum2(min1)
        self.pick_count += 1

    def on_btnDisconnectivity_graph_clicked(self, clicked=None):
        """show the disconnectivity graph 
        
        it is interactive, so that when you click on an end point
        that minima is selected
        """
        if clicked is None: return

        if not hasattr(self, "dgraph_dlg"):
            self.dgraph_dlg = DGraphDialog(self.system.database, parent=self)
            self.dgraph_dlg.dgraph_widget.minimum_selected.connect(self.on_minimum_picked)
        self.dgraph_dlg.rebuild_disconnectivity_graph()
        self.dgraph_dlg.show()

    
    def on_btnShowGraph_clicked(self, clicked=None):
        """ show the graph of minima and transition states 
        
        make it interactive, so that when you click on a point
        that minima is selected
        """
        if clicked is None: return
        self.pick_count = 0
        from pygmin.gui.graph_viewer import GraphViewDialog
        if not hasattr(self, "graphview"):
            self.graphview = GraphViewDialog(self.system.database, parent=self, app=self.app)
            self.graphview.widget.on_minima_picked.connect(self.on_minimum_picked)
        self.graphview.show()
        self.graphview.widget.make_graph()
        self.graphview.widget.show_graph()
        return
        
    def on_pushNormalmodesMin_clicked(self, clicked=None):
        if clicked is None: return
        if not hasattr(self, "normalmode_explorer"):
            self.normalmode_explorer = NormalmodeBrowser(self, self.system, self.app)
        min1 = self.ui.ogl_main.minima[1]
        self.normalmode_explorer.set_coords(min1.coords)
        self.normalmode_explorer.show()
    
    def get_selected_ts(self):
        ts = self.ts_selected
        if ts is None:
            raise Exception("you must select a transition state first")
        return ts
    
    def on_pushNormalmodesTS_clicked(self, clicked=None):
        if clicked is None: return
        if not hasattr(self, "normalmode_explorer"):
            self.normalmode_explorer = NormalmodeBrowser(self, self.system, self.app)
        ts = self.get_selected_ts()
        self.normalmode_explorer.set_coords(ts.coords)
        self.normalmode_explorer.show()
    
    def _sort_ts(self):
        self.need_sorting_ts = True
        self._sort_lists()
    
    def _sort_minima(self):
        self.need_sorting_minima = True
        self._sort_lists()
    
    def _sort_lists(self):
        """calling this function indicates that the lists need sorting
        
        since sorting can be a *huge* bottleneck for large lists we try to do it
        as infrequently as possible.  Here we wait several seconds too see
        if more minima are added before sorting the lists
        """
        self.need_sorting = True
        if not self._sort_timer.isActive():
            self._sort_timer.start(2000)

    def _delayed_sort(self):
        if not self.need_sorting:
            self.need_sorting_minima = False
            self.need_sorting_ts = False
            print "sorting not needed"
            return
        try:
            # the system params flag flag gui._sort_lists can optionally
            # be set to False to indicate not to sort the lists.  This is
            # useful if many minima are added at once, e.g. at the end of
            # a connect run.  The flag should be set to True when all the
            # minima are added. 
            s = self.system.params.gui._sort_lists
#            print "self.system.params.gui._sort_lists", s
            if not s:
                # wait a few seconds and call this function again
                if not self._sort_timer.isActive():
                    self._sort_timer.start(2000)
                print "delaying sort"
                return
        except AttributeError:
            pass
        
#        print "sorting lists"
        if self.need_sorting_minima:
            self.minima_list_model.sort(0)
        if self.need_sorting_ts:
            self.ts_list_model.sort(0)
        self.need_sorting = False
        self.need_sorting_minima = False
        self.need_sorting_ts = False
        self._sort_timer.stop()
#        print "done sorting lists"
        
                    
    
    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
            
        self.minima_list_model.appendRow(MinimumStandardItem(minimum))
        if sort_items:
            self._sort_minima()
                
    def RemoveMinimum(self, minimum):
        """remove a minimum from self.minima_list_model"""
        minid = minimum._id
        items = self.minima_list_model.findItems('*', QtCore.Qt.MatchWildcard)
        for i in items:
#            print "item", i.minimum._id, minid, minimum._id
            if i.minimum._id == minid:
#                print "taking item", i.minimum._id
                self.minima_list_model.takeRow(i.row())
                break

    def NewTS(self, ts, sort=True):
        """add new transition state, or list of transition states"""
        try:
            len(ts)
            is_iterable = True
        except TypeError:
            is_iterable = False
        if is_iterable:   
            tslist = ts
        else: 
            tslist = [ts]
        for ts in tslist:
            tsitem = TransitionStateStandardItem(ts)
            self.ts_list_model.appendRow(tsitem)
        if sort:
            self._sort_ts()

    def RemoveTS(self, ts):
        """remove transition state"""
        raise Exception("removing transition states not implemented yet")
        obj = self.ui.list_TS
        tsid = id(ts)
        itms = self.ui.list_TS.findItems('*', QtCore.Qt.MatchWildcard)
        for i in itms:
            if i.tsid == tsid:
                obj.takeItem(obj.row(i))

  
                     
    def on_btn_start_basinhopping_clicked(self, clicked=None):
        if clicked is None: return
        db = self.system.database
        self.system.database = None
        self.bhrunner = BHRunner(self.system)
        self.bhrunner.start()
        self.system.database = db
        
#    def tsSearch(self):
#        ts = self.system.findTS(self.ui.oglTS.coords[1])
#        self.transition = [ts[1][0], ts[0][0], ts[2][0]]

#    def showFrameTS(self, i):
#        if self.transition:
#            self.ui.oglTS.setCoords(self.transition[i])

    def on_action_delete_minimum_triggered(self, checked=None):
        if checked is None: return
        min1 = self.ui.ogl_main.minima[1]
        ret = QtGui.QMessageBox.question(self, "Deleting minima", 
                                   "Do you want to delete minima %d with energy %g"%(min1._id, min1.energy), 
                                   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel)
        if(ret == QtGui.QMessageBox.Ok):
            print "deleting minima"
            print "deleting minimum", min1._id, min1.energy
            self.RemoveMinimum(min1)
            self.system.database.removeMinimum(min1)
    
    def on_btnConnect_clicked(self, clicked=None):
        if clicked is None: return
        return self._doubleEndedConnect(reconnect=False)

    def on_btnReconnect_clicked(self, clicked=None):
        if clicked is None: return
        return self._doubleEndedConnect(reconnect=True)

    def _doubleEndedConnect(self, reconnect=False, min1min2=None):
        """
        launch a double ended connect run to connect the two selected minima.

        If the minima are not connected, or reconnect is True, launch a connect browser
        in  a separate window.  Else just show the path in the OGL viewer
        """
        # determine which minima to connect
        if min1min2 is None:
            min1, min2 = self.get_selected_minima()
        else:
            min1, min2 = min1min2
        database = self.system.database
        if not reconnect:
            # check if the minima are already connected
            double_ended_connect = self.system.get_double_ended_connect(min1, min2, database, 
                                                   fresh_connect=False, verbosity=0)
            if double_ended_connect.graph.areConnected(min1, min2):
                print "minima are already connected.  loading smoothed path in viewer"
                mints, S, energies = double_ended_connect.returnPath()
                clist = [m.coords for m in mints]
                smoothpath = self.system.smooth_path(clist)
                
                coords = np.array(smoothpath)
                self.nebcoords = coords
                self.nebenergies = np.array(energies)
                print "setting path in oglPath"
                self.ui.oglPath.setCoordsPath(coords)#, labels)
#                self.ui.oglPath.setCoords(coords[0,:], 1)
#                self.ui.oglPath.setCoords(None, 2)
#                self.ui.sliderFrame.setRange(0, coords.shape[0]-1)
                
                if self.usepymol:
                    self.pymolviewer.update_coords(self.nebcoords, index=1, delete_all=True)
                
                return

                
        # make the connect viewer
        
        
        decviewer = ConnectViewer(self.system, self.system.database, min1, min2, parent=self, app=self.app)
        
        print "starting double ended"
        decviewer.show()
        decviewer.start()
        
        # store pointers
        self.double_ended_connect_runs.append(decviewer)
     
    def on_btn_connect_in_optim_clicked(self, clicked=None):
        """spawn an OPTIM job and retrieve the minima and transition states 
        it finds"""
        if clicked is None: return
        min1, min2 = self.get_selected_minima()
#        existing_minima = set(self.system.database.minima())
        spawner = self.system.get_optim_spawner(min1.coords, min2.coords)
        spawner.run()
        db = self.system.database
        newminima, newts = spawner.load_results(self.system.database)
#        for m in newminima:
#            if m not in existing_minima:
#                self.NewMinimum(m)
            
        #now use DoubleEndedConnect to test if they are connected
        graph = Graph(db)
        if graph.areConnected(min1, min2):
            #use double ended connect to draw the interpolated path
            #this is ugly
            self._doubleEndedConnect(reconnect=False, min1min2=(min1, min2))
                                                       

    def _merge_minima(self, min1, min2):
        mindist = self.system.get_mindist()
        dist, x1, x2 = mindist(min1.coords, min2.coords)
        query  = "Do you want to merge minimum %d with energy %g" %(min1._id, min1.energy)
        query += "                with minimum %d with energy %g" %(min2._id, min2.energy)
        query += "    separated by distance %g" % (dist)
        ret = QtGui.QMessageBox.question(self, "Merging minima", 
                                   query, 
                                   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel)
        if(ret == QtGui.QMessageBox.Ok):
            m1, m2 = min1, min2
            if m1._id > m2._id:
                m1, m2 = m2, m1
            print "merging minima", m1._id, m2._id#, ": minimum", m2._id, "will be deleted"
            self.system.database.mergeMinima(m1, m2)
            self.RemoveMinimum(m2)

    def on_action_merge_minima_triggered(self, checked=None):
        if checked is None: return
        min1, min2 = self.get_selected_minima()
        self._merge_minima(min1, min2)



#    def launch_connect_explorer(self):
#        coords1, coords2 = self.get_selected_coords()
#
#        if not hasattr(self, "local_connect_explorer"):
#            self.local_connect_explorer = ConnectExplorerDialog(self.system)
#            self.local_connect_explorer.nebwgt.process_events.connect(self.processEvents)
#        self.local_connect_explorer.show()
#        self.local_connect_explorer.createNEB(coords1, coords2)
#        self.local_connect_explorer.runNEB()


    def on_btn_close_all_clicked(self, checked=None):
        if checked is None: return
        print "closing all windows"
        for dv in self.double_ended_connect_runs:
            dv.hide()
#            del dv
        self.double_ended_connect_runs = []

        try:
            self.local_connect_explorer.hide()
            del self.local_connect_explorer
        except AttributeError: pass

        try:        
            self.dgraph_dlg.hide()
            del self.dgraph_dlg
        except AttributeError: pass
        
        try:
            self.nebexplorer.hide()
            del self.nebexplorer
        except AttributeError: pass

    def on_btn_connect_all_clicked(self, checked=None):
        if checked is None: return
        from pygmin.gui.connect_all import ConnectAllDialog
#        if hasattr(self, "connect_all"):
#            if not self.connect_all.isVisible():
#                self.connect_all.show()
#            if not self.connect_all.is_running()
        self.connect_all = ConnectAllDialog(self.system, self.system.database, 
                                            parent=self, app=self.app)
        self.connect_all.show()
        self.connect_all.start()
        
    def on_pushTakestepExplorer_clicked(self):
        if not hasattr(self, "takestep_explorer"):
            self.takestep_explorer = TakestepExplorer(parent=self, system = self.system, app = self.app,
                                                      database = self.system.database)
            
        self.takestep_explorer.show()
        
#def refresh_pl():
    #pl.pause(0.000001)    
    
def run_gui(system, db=None):
    """
    The top level function that will launch the gui for a given system
    
    Parameters
    ----------
    system : System class
        A pygmin system, derived from BaseSystem.  All information 
        about the system is in this class.
    db : str, optional
        connect to the database at this file location
        
    """
    app = QtGui.QApplication(sys.argv)
    
    sys.excepthook = excepthook

    myapp = MainGUI(app, system)
    if db is not None:
        myapp.connect_db(db)
        
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
    myapp.show()
    sys.exit(app.exec_()) 
       
#def run_gui(systemtype):
#    app = QtGui.QApplication(sys.argv)
#    import pylab as pl
#    myapp = MainGUI(systemtype)
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
#    
#    myapp.show()
#    sys.exit(app.exec_())
