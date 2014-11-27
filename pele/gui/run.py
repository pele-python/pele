import matplotlib
matplotlib.use("QT4Agg")
import traceback    
import sys
import numpy as np

from PyQt4 import QtCore, QtGui

from pele.landscape import TSGraph
from pele.storage import Database
from pele.utils.events import Signal
from pele.config import config

from pele.gui.MainWindow import Ui_MainWindow 
from pele.gui.bhrunner import BHManager
from pele.gui.dlg_params import DlgParams
from pele.gui.ui.dgraph_dlg import DGraphDialog
from pele.gui.connect_run_dlg import ConnectViewer
from pele.gui.takestep_explorer import TakestepExplorer
from pele.gui.normalmode_browser import NormalmodeBrowser
from pele.gui._list_views import ListViewManager
from pele.gui._cv_viewer import HeatCapacityViewer
from pele.gui._rate_gui import RateViewer
from pele.gui.graph_viewer import GraphViewDialog



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
        raise ex_value

class MySelection(object):
    """keep track of which minima have been selected and whether those coordinates have been modified
    """
    def __init__(self):
        self.minimum1 = None
        self.minimum2 = None
        self.coords1 = None
        self.coords2 = None


class MainGUI(QtGui.QMainWindow):
    """
    this is the main class for the pele gui
    
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

        self.minima_selection = MySelection()
        self.on_minimum_1_selected = Signal()
        self.on_minimum_2_selected = Signal()
 
        # set up the list manager
        self.list_manager = ListViewManager(self)
        
        # define the system
        self.NewSystem()

        # finish setting up the list manager (this must be done after NewSystem() is called)
        self.list_manager.finish_setup()

        # try to load the pymol viewer.
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
            
        if not self.usepymol:
            # note: glutInit() must be called exactly once.  pymol calls it
            # during pymol.finish_launching(), so if we call it again it will
            # give an error. On the other hand, if we're not using pymol we 
            # must call it.
            from OpenGL.GLUT import glutInit
            glutInit()
        
        self.bhmanager = None
        
    def NewSystem(self):
        """
        this is called to initialize the system with a database
        """
        self.system = self.systemtype()
        self.connect_db()
            
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
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Select database', '.')
        if len(filename) > 0:
            self.connect_db(filename)

    def connect_db(self, database=":memory:"):
        """
        connect to an existing database at location filename
        """
        self.list_manager.clear()

        # note: database can be either Database, or string, or QString        
        if isinstance(database, Database):
            self.system.database = database
        else:
            self.system.database = self.system.create_database(db=database)
        # add minima to listWidged.  do sorting after all minima are added
        for minimum in self.system.database.minima():
            self.NewMinimum(minimum, sort_items=False)
        self.list_manager._sort_minima()
        self.NewTS(self.system.database.transition_states(order_energy=True))
        self.list_manager.resize_columns_minima()
        self.list_manager.resize_columns_ts()

        self.system.database.on_minimum_added.connect(self.NewMinimum)
        self.system.database.on_minimum_removed(self.RemoveMinimum)
        self.system.database.on_ts_added.connect(self.NewTS)
        self.system.database.on_ts_removed.connect(self.RemoveTS)
    


    def SelectMinimum(self, minimum, set_selected=True):
        """when you click on a minimum in the basinhopping tab
        """
#        print "selecting minimum", minimum._id, minimum.energy
        if set_selected:
            self.list_manager._select_main(minimum)
            return
        self.ui.ogl_main.setSystem(self.system)
        self.ui.ogl_main.setCoords(minimum.coords)
        self.ui.ogl_main.setMinimum(minimum)
        self.ui.oglTS.setSystem(self.system)
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=1, delete_all=True)

    def _SelectMinimum1(self, minimum, set_selected=True):
        """set the first minimum displayed in the connect tab"""
        if set_selected:
            self.list_manager._select1(minimum)
            return
        print "selecting minimum 1:", minimum._id, minimum.energy
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=1)
#        self.ui.oglPath.setMinimum(minimum, index=1)
        self.minima_selection.minimum1 = minimum
        self.minima_selection.coords1 = minimum.coords
        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=1)
        
        self.on_minimum_1_selected(minimum)
        
    def _SelectMinimum2(self, minimum, set_selected=True):
        """set the second minimum displayed in the connect tab"""
        if set_selected:
            self.list_manager._select2(minimum)
            return
        print "selecting minimum 2:", minimum._id, minimum.energy
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=2)
#        self.ui.oglPath.setMinimum(minimum, index=2)
        self.minima_selection.minimum2 = minimum
        self.minima_selection.coords2 = minimum.coords

        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=2)
        
        self.on_minimum_2_selected(minimum)
    
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
        self.nebexplorer.new_neb(coords1, coords2, run=False)
        
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
            self.dgraph_dlg.dgraph_widget.minimum_selected.connect(self.SelectMinimum)
        self.dgraph_dlg.rebuild_disconnectivity_graph()
        self.dgraph_dlg.show()

    
    def on_btnShowGraph_clicked(self, clicked=None):
        """ show the graph of minima and transition states 
        
        make it interactive, so that when you click on a point
        that minima is selected
        """
        if clicked is None: return
        self.pick_count = 0
        if not hasattr(self, "graphview"):
            self.graphview = GraphViewDialog(self.system.database, parent=self, app=self.app)
            self.graphview.widget.on_minima_picked.connect(self.on_minimum_picked)
            self.graphview.widget.on_minima_picked.connect(self.SelectMinimum)
        self.graphview.show()
        self.graphview.widget.make_graph()
        try:
            m1, m2 = self.get_selected_minima()
            self.graphview.widget._show_minimum_energy_path(m1, m2)
        except Exception:
            self.graphview.widget.show_graph()
        
    def on_pushNormalmodesMin_clicked(self, clicked=None):
        if clicked is None: return
        if not hasattr(self, "normalmode_explorer"):
            self.normalmode_explorer = NormalmodeBrowser(self, self.system, self.app)
        min1 = self.ui.ogl_main.minima[1]
        if min1 is None:
            raise RuntimeError("you must select a minimum first")
        self.normalmode_explorer.set_coords(min1.coords)
        self.normalmode_explorer.show()
    
    def on_pushNormalmodesTS_clicked(self, clicked=None):
        if clicked is None: return
        if not hasattr(self, "normalmode_explorer"):
            self.normalmode_explorer = NormalmodeBrowser(self, self.system, self.app)
        ts = self.list_manager.get_selected_ts()
        if ts is None:
            raise RuntimeError("you must select a transition state first")
        self.normalmode_explorer.set_coords(ts.coords)
        self.normalmode_explorer.show()

    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
        self.list_manager.NewMinimum(minimum, sort_items=sort_items)
                
    def RemoveMinimum(self, minimum):
        """remove a minimum from self.minima_list_model"""
        self.list_manager.RemoveMinimum(minimum)

    def NewTS(self, ts, sort=True):
        """add new transition state, or list of transition states"""
        self.list_manager.NewTS(ts, sort=sort)

    def RemoveTS(self, ts):
        """remove transition state"""
        raise Exception("removing transition states not implemented yet")
        obj = self.ui.list_TS
        tsid = id(ts)
        itms = self.ui.list_TS.findItems('*', QtCore.Qt.MatchWildcard)
        for i in itms:
            if i.tsid == tsid:
                obj.takeItem(obj.row(i))

    def set_basinhopping_number_alive(self, nalive):
        """set the label that shows how many basinhopping processes are alive"""
        self.ui.label_bh_nproc.setText("%d B.H. processes" % nalive)
  
    def on_btn_start_basinhopping_clicked(self, clicked=None):
        """this is run when the start basinhopping button is clicked"""
        if clicked is None: return
        # set up the basinhopping manager if not already done
        if self.bhmanager is None:
            self.bhmanager = BHManager(self.system, self.system.database,
                                       on_number_alive_changed=self.set_basinhopping_number_alive)
        # get the number of steps from the input box
        nstepsstr = self.ui.lineEdit_bh_nsteps.text()
        nsteps = None
        try:
            nsteps = int(nstepsstr)
        except ValueError:
            # ignore the text if it is the default text
            if "steps" not in nstepsstr:
                sys.stderr.write("can't convert %s to integer\n" % nstepsstr)
        
        # start a basinhopping run
        self.bhmanager.start_worker(nsteps=nsteps)

    def on_btn_stop_basinhopping_clicked(self, clicked=None):
        if clicked is None: return
        self.bhmanager.kill_all_workers()

    def on_action_delete_minimum_triggered(self, checked=None):
        if checked is None: return
        min1 = self.ui.ogl_main.minima[1]
        ret = QtGui.QMessageBox.question(self, "Deleting minima", 
                                   "Do you want to delete minima %d with energy %g"%(min1._id, min1.energy), 
                                   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel)
        if ret == QtGui.QMessageBox.Ok:
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
                self.ui.oglPath.setCoordsPath(coords)
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
            
        # now use DoubleEndedConnect to test if they are connected
        graph = TSGraph(db)
        if graph.areConnected(min1, min2):
            # use double ended connect to draw the interpolated path
            # this is ugly
            self._doubleEndedConnect(reconnect=False, min1min2=(min1, min2))
                                                       

    def _merge_minima(self, min1, min2):
        mindist = self.system.get_mindist()
        dist, x1, x2 = mindist(min1.coords, min2.coords)
        query  = "Do you want to merge minimum %d with energy %g" %(min1._id, min1.energy)
        query += "                with minimum %d with energy %g" %(min2._id, min2.energy)
        query += "    separated by distance %g" % dist
        ret = QtGui.QMessageBox.question(self, "Merging minima", 
                                   query, 
                                   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel)
        if ret == QtGui.QMessageBox.Ok:
            m1, m2 = min1, min2
            if m1._id > m2._id:
                m1, m2 = m2, m1
            print "merging minima", m1._id, m2._id
            self.system.database.mergeMinima(m1, m2)
            self.RemoveMinimum(m2)

    def on_action_merge_minima_triggered(self, checked=None):
        if checked is None: return
        min1, min2 = self.get_selected_minima()
        self._merge_minima(min1, min2)

    def on_action_compute_thermodynamic_info_triggered(self, checked=None):
        if checked is None: return
        def on_done(): print "done computing thermodynamic info"
        self._on_done = on_done # because on_finish stores a weak reference
        self.compute_thermodynamic_information(on_finish=self._on_done )

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
        
        try:
            self.rate_viewer.hide()
            del self.rate_viewer
        except AttributeError: pass

    def on_btn_connect_all_clicked(self, checked=None):
        if checked is None: return
        from pele.gui.connect_all import ConnectAllDialog
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
    
    def on_btn_heat_capacity_clicked(self, clicked=None):
        if clicked is None: return
        self.cv_viewer = HeatCapacityViewer(self.system, self.system.database, parent=self)
        self.cv_viewer.show()
        self.cv_viewer.rebuild_cv_plot()
    
    def compute_thermodynamic_information(self, on_finish=None):
        """compute thermodynamic information for minima and ts in the background
        
        call on_finish when the calculation is done
        """
        # TODO: deal carefuly with what will happen if this is called again
        # before the first calculation is done.  if self.thermo_worker is overwritten will
        # the first calculation stop?
        from pele.gui._cv_viewer import GetThermodynamicInfoParallelQT
        self.thermo_worker = GetThermodynamicInfoParallelQT(self.system, self.system.database, npar=1)
        if on_finish is not None:
            self.thermo_worker.on_finish.connect(on_finish)
        self.thermo_worker.start()
        njobs = self.thermo_worker.njobs
        print "calculating thermodynamics for", njobs, "minima and transition states" 
        

#    def _compute_rates(self, min1, min2, T=1.):
#        """compute rates without first calculating thermodynamics
#        """
#        print "computing rates at temperature T =", T
#        tslist = [ts for ts in self.system.database.transition_states() 
#                  if ts.fvib is not None]
#        rcalc = RateCalculation(tslist, [min1], [min2], T=T)
#        r12, r21 = rcalc.compute_rates()
#        print "rate from", min1._id, "to", min2._id, "=", r12
#        print "rate from", min2._id, "to", min1._id, "=", r21
#
#    def compute_rates(self, min1, min2, T=1.):
#        """compute the transition rate from min1 to min2 and vice versa"""
#        def on_finish():
#            print "thermodynamic calculation finished"
#            self._compute_rates(min1, min2)
#        self._on_finish_thermo_reference = on_finish # so it doeesn't get garbage collected
#        self.compute_thermodynamic_information(on_finish=on_finish)
    
    def on_btn_rates_clicked(self, clicked=None):
        if clicked is None: return
        if not hasattr(self, "rate_viewer"):
            m1, m2 = self.minima_selection.minimum1, self.minima_selection.minimum2
            self.rate_viewer = RateViewer(self.system, self.system.database, parent=self)
            if m1 is not None:
                self.rate_viewer.update_A(m1)
            if m2 is not None:
                self.rate_viewer.update_B(m2)
            self.on_minimum_1_selected.connect(self.rate_viewer.update_A)
            self.on_minimum_2_selected.connect(self.rate_viewer.update_B)
        self.rate_viewer.show()

        
def run_gui(system, db=None, application=None):
    """
    The top level function that will launch the gui for a given system
    
    Parameters
    ----------
    system : System class
        A pele system, derived from BaseSystem.  All information 
        about the system is in this class.
    db : pele database or string, optional
        connect to this database or the database at this file location
    application : QApplication
        Use this QApplication object rather than creating a new one
        
    """
    if application is None:
        application = QtGui.QApplication(sys.argv)
    
    sys.excepthook = excepthook

    myapp = MainGUI(application, system)
    if db is not None:
        myapp.connect_db(db)
        
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
    myapp.show()
    sys.exit(application.exec_()) 
       
