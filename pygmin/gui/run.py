import matplotlib
matplotlib.use("QT4Agg")
    
import pylab as pl    
from PyQt4 import QtCore, QtGui
import MainWindow 
from collections import deque
import sys
import bhrunner
import copy
import numpy as np
#import matplotlib.pyplot as plt

from pygmin.storage import Database
from pygmin.landscape import Graph
from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from dlg_params import DlgParams
from pygmin.config import config

global pick_count

def no_event(*args, **kwargs):
    return

class pick_event(object):
    """
    define a matplotlib pick_event event
    """
    def __init__(self, artist):
        pass

class QMinimumInList(QtGui.QListWidgetItem):
    
    def __init__(self, minimum):
        text="%.4f (%d)"%(minimum.energy, minimum._id)
        QtGui.QListWidgetItem.__init__(self, text)
        
    def setCoords(self, coords):
        self.coords = coords
    def setMinimum(self, minimum):
        self.minid = id(minimum)
        self.minimum = minimum
    def __lt__(self, item2):
        #sort the energies in the list lowest to highest
        return self.minimum.energy > item2.minimum.energy

class MyForm(QtGui.QMainWindow):
    def __init__(self, app, systemtype, parent=None):
        QtGui.QWidget.__init__(self)
        self.ui = MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.listMinima = [
                           self.ui.listWidget,
                           self.ui.listMinima1,
                           self.ui.listMinima2,
                           self.ui.listFrom
                           ]
        self.systemtype = systemtype
        self.NewSystem()
        self.transition=None
        self.app = app
        
        #try to load the pymol viewer.  
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
        self.system = self.systemtype()
        db = self.system.create_database()
        self.system.database = db
        self.system.database.on_minimum_added.connect(self.NewMinimum)
        self.system.database.on_minimum_removed.connect(self.RemoveMinimum)
        for l in self.listMinima:
            l.clear()
            
    def edit_params(self):
        self.paramsdlg = DlgParams(self.system.params)
        self.paramsdlg.show()
        
    #def save(self):
    #    import pickle
    #    filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '.')
    #    output = open(filename, "w")
    #    pickle.dump(self.system.storage, output)
       
    def connect(self):
        """
        connect to an existing database
        """
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Open File', '.')
        self.connect_db(filename)

    def connect_db(self, filename):
        db = self.system.create_database(db=filename)
        self.system.database = db
        #add minima to listWidged.  do sorting after all minima are added
        for minimum in self.system.database.minima():
            self.NewMinimum(minimum, sort_items=False)
        for obj in self.listMinima:
            obj.sortItems(1)

        self.system.database.on_minimum_added.connect(self.NewMinimum)
        self.system.database.on_minimum_removed(self.RemoveMinimum)
        
    def SelectMinimum(self, item):
        """when you click on a minimum in the basinhopping tab"""
        print "selecting minimum", item.minimum._id, item.minimum.energy
        self.ui.widget.setSystem(self.system)
        self.ui.widget.setCoords(item.coords)
        self.ui.widget.setMinimum(item.minimum)
        self.ui.oglTS.setSystem(self.system)
        self.ui.oglTS.setCoords(item.coords)
        if self.usepymol:
            self.pymolviewer.update_coords([item.coords], index=1, delete_all=True)
        
    def _SelectMinimum1(self, minimum):
        """by minimum"""
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=1)
        self.ui.oglPath.setMinimum(minimum, index=1)
        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=1)


    def SelectMinimum1(self, item):
        """called by the ui"""
        print "selecting minimum 1", item.minimum._id, item.minimum.energy
        return self._SelectMinimum1(item.minimum)
 
    def _SelectMinimum2(self, minimum):
        """by minimum"""
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(minimum.coords, index=2)
        self.ui.oglPath.setMinimum(minimum, index=2)
        self.neb = None
        if self.usepymol:
            self.pymolviewer.update_coords([minimum.coords], index=2)


    def SelectMinimum2(self, item):
        """called by the ui"""
        print "selecting minimum 2", item.minimum._id, item.minimum.energy
        return self._SelectMinimum2(item.minimum)
    
    
    def Invert(self):
        """invert the coordinates"""
        coords2 = self.ui.oglPath.coords[2]
        self.ui.oglPath.setCoords(-coords2, 2)
        if self.usepymol:
            self.pymolviewer.update_coords([-coords2], index=2)

    
    def AlignMinima(self):
        """use mindist to align the minima"""
        coords1 = self.ui.oglPath.coords[1]
        coords2 = self.ui.oglPath.coords[2]
        align = self.system.get_mindist()
        pot = self.system.get_potential()
        print "energy before alignment", pot.getEnergy(coords1), pot.getEnergy(coords2)
        dist, coords1, coords2 = align(coords1, coords2)
        print "energy after alignment", pot.getEnergy(coords1), pot.getEnergy(coords2)
        self.ui.oglPath.setCoords(coords1, 1)
        self.ui.oglPath.setCoords(coords2, 2)
        print "best alignment distance", dist
        if self.usepymol:
            self.pymolviewer.update_coords([coords1], index=1)
            self.pymolviewer.update_coords([coords2], index=2)
    
    def ConnectMinima(self):
        """do an NEB run (not a connect run).  Don't find best alignment first"""
        coords1 = self.ui.oglPath.coords[1]
        coords2 = self.ui.oglPath.coords[2]
        min1 = self.ui.oglPath.minima[1]
        min2 = self.ui.oglPath.minima[2]
        
        #use the functions in DoubleEndedConnect to set up the NEB in the proper way
        double_ended = self.system.get_double_ended_connect(min1, min2, 
                                                            self.system.database, 
                                                            fresh_connect=True)
        
        # the following lines implement the ability to follow the status of the
        # NEB in real time.  the callback is passed to the NEB and plots the
        # last several values of the energies of the path.  The implementation
        # is quite simple and could easily be a lot faster.
        # set follow_neb=False to turn this off
        class NEBCallback(object):
            def __init__(self, app, frq=30, nplots=3):
                self.count = 0
                self.nplots = nplots
                self.data = deque()
                self.frq = frq
                self.app = app
                
            def __call__(self, energies=None, coords=None, **kwargs):
                self.count += 1
                if self.count % self.frq == 1:
                    self.data.append(energies.copy())
                    if len(self.data) > self.nplots:
                        self.data.popleft()
                    pl.clf()
                    for E in self.data:
                        line, = pl.plot(E, "o-")
                        # note: if we save the line and use line.set_ydata(E)
                        # this would be a lot faster, but we would have to keep
                        # track of the y-axis limits manually
                    pl.draw()
                    self.app.processEvents()

        follow_neb = True
        if follow_neb:
            pl.ion()
            pl.clf()
            neb_callback = NEBCallback(self.app)
        else:
            neb_callback = no_event
        
        local_connect = double_ended._getLocalConnectObject()
        self.neb =  local_connect._getNEB(self.system.get_potential(),
                                          coords1, coords2, event=neb_callback,
                                          verbose=True,
                                          **local_connect.NEBparams)        
        
        self.neb.optimize()
        self.nebcoords = self.neb.coords
        self.nebenergies = self.neb.energies
        self.ui.oglPath.setCoords(self.neb.coords[0,:], 1)
        self.ui.oglPath.setCoords(None, 2)
        self.ui.sliderFrame.setRange(0,self.neb.coords.shape[0]-1)
        if self.usepymol:
            self.pymolviewer.update_coords(self.nebcoords, index=1, delete_all=True)

    def showFrame(self, i):
        if hasattr(self, "nebcoords"):
            self.ui.oglPath.setCoords(self.nebcoords[i,:])
    
    def show_disconnectivity_graph(self):
        """show the disconnectivity graph 
        
        make it interactive, so that when you click on an end point
        that minima is selected
        """
        import pylab as pl
        pl.ion()
        pl.clf()
        ax = pl.gca()
        fig = pl.gcf()
        
        ax.grid(True)
        

        graphwrapper = Graph(self.system.database)
        dg = DisconnectivityGraph(graphwrapper.graph, subgraph_size=2)
        dg.calculate()
        
        #draw minima as points
        xpos, minima = dg.get_minima_layout()
        energies = [m.energy for m in minima]
        points = ax.scatter(xpos, energies, picker=5)
        
        #draw line segments connecting minima
        line_segments = dg.line_segments
        for x, y in line_segments:
            ax.plot(x, y, 'k')
        
        
        #define what happens when a point is clicked on
        global pick_count
        pick_count = 0
        def on_pick(event):
            if event.artist != points:
#                print "you clicked on something other than a node"
                return True
            thispoint = event.artist
            ind = event.ind[0]
            min1 = minima[ind]
            print "you clicked on minimum with id", min1._id, "and energy", min1.energy
            global pick_count
            #print pick_count
            pick_count += 1
            if (pick_count % 2) == 0:
                self._SelectMinimum1(min1)
            else:
                self._SelectMinimum2(min1)
        fig = pl.gcf()
        cid = fig.canvas.mpl_connect('pick_event', on_pick)

        pl.show()
        

    
    def show_graph(self):
        """ show the graph of minima and transition states 
        
        make it interactive, so that when you click on a point
        that minima is selected
        """
        import pylab as pl
        import networkx as nx
        pl.ion()
        pl.clf()
        ax = pl.gca()
        fig = pl.gcf()
        
        #get the graph object, eliminate nodes without edges
        graphwrapper = Graph(self.system.database)
        graph = graphwrapper.graph
        degree = graph.degree()
        nodes = [n for n, nedges in degree.items() if nedges > 0]
        graph = graph.subgraph(nodes)
        
        #get the layout of the nodes from networkx
        layout = nx.spring_layout(graph)
        layoutlist = layout.items()
        xypos = np.array([xy for n, xy in layoutlist])
        #color the nodes by energy
        e = np.array([m.energy for m, xy in layoutlist])
        #plot the nodes
        points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                            s=8**2, c=e, cmap=pl.cm.autumn)
        fig.colorbar(points)
        #label the nodes
        ids = [n._id for n, xy in layoutlist]
        for i in range(len(ids)):
            ax.annotate( ids[i], xypos[i] )
        
        
    
        #plot the edges as lines
        for u, v in graph.edges():
            line = np.array([layout[u], layout[v]])
            ax.plot(line[:,0], line[:,1], '-k')
        
        #scale the axes so the points are not cutoff
        xmin, ymin = np.min(xypos, 0)
        xmax, ymax = np.max(xypos, 0)
        dx = (xmax - xmin)*.1
        dy = (ymax - ymin)*.1
        ax.set_xlim([xmin-dx, xmax+dx])
        ax.set_ylim([ymin-dy, ymax+dy])
        
        global pick_count
        pick_count = 0

        def on_pick(event):
            if event.artist != points:
#                print "you clicked on something other than a node"
                return True
            thispoint = event.artist
            ind = event.ind[0]
            min1 = layoutlist[ind][0]
            print "you clicked on minimum with id", min1._id, "and energy", min1.energy
            global pick_count
            #print pick_count
            pick_count += 1
            if (pick_count % 2) == 0:
                self._SelectMinimum1(min1)
            else:
                self._SelectMinimum2(min1)

        fig = pl.gcf()
        cid = fig.canvas.mpl_connect('pick_event', on_pick)
        
        #ids = dict([(n, n._id) for n in nodes])
        #nx.draw(graph, labels=ids, nodelist=nodes)
        #ax.draw()
#        pl.draw()
        pl.show()
        
    
    def showEnergies(self):
        """plot the energies from NEB or connect
        
        don't clear the previous plot so we can overlay multiple plots
        """
        #note: this breaks if pylab isn't a local import.  I don't know why
        import pylab as pl
        pl.ion()
        ax = pl.gca()
        ax.plot(self.nebenergies, "-")
        points = ax.scatter(range(len(self.nebenergies)), self.nebenergies, picker=5, label="energies")
        
        if False:
            #start to make the energies interactive.  I would like it to be such that
            #when you click on a point, that structure gets selected. But the structures
            #in the NEB are not Minimum objects, so they can't be selected using the current
            #setup
            def on_pick(event):
                if event.artist != points:
                    print "you clicked on something other than a node"
                    return True
                ind = event.ind[0]
            #    min1 = layoutlist[ind][0]
                yvalue = self.nebenergies[ind]
                print "you clicked on a configuration with energy", yvalue

            
            fig = pl.gcf()
            cid = fig.canvas.mpl_connect('pick_event', on_pick)


        
        if False: #show climbing images
            neb = self.neb
            cl=[]
            en=[]
            for i in xrange(len(neb.energies)):
                if(neb.isclimbing[i]):
                    print "climbing image :", i, neb.energies[i]
                    cl.append(i)
                    en.append(neb.energies[i])
                    
            pl.plot(cl, en, "s", label="climbing images", markersize=10, markerfacecolor="none", markeredgewidth=2)
        
        pl.legend(loc='best')
        pl.show()
     
    def NewMinimum(self, minimum, sort_items=True):
        """ add a new minimum to the system """
        E=minimum.energy
        minid=id(minimum)
        coords=minimum.coords
        for obj in self.listMinima:
            item = QMinimumInList(minimum)
            item.setCoords(coords)
            item.setMinimum(minimum)
            obj.addItem(item)    
            if sort_items:
                obj.sortItems(1)
                
    def RemoveMinimum(self, minimum):
        minid = id(minimum)
        for obj in self.listMinima:
            itms = obj.findItems('*', QtCore.Qt.MatchWildcard)
            for i in itms:
                if(i.minid == minid):
                    obj.takeItem(obj.row(i))
                        
    def StartBasinHopping(self):
        db = self.system.database
        self.system.database = None
        self.bhrunner = bhrunner.BHRunner(self.system)
        self.bhrunner.start()
        self.system.database = db
        
    def tsSearch(self):
        import numpy as np
        ts = self.system.findTS(self.ui.oglTS.coords[1])
        self.transition = [ts[1][0], ts[0][0], ts[2][0]]

    def showFrameTS(self, i):
        if(self.transition):
            self.ui.oglTS.setCoords(self.transition[i])

    def selectTransition(self):
        pass

    def delete_minimum(self):
        min1 = self.ui.widget.minima[1]
        ret = QtGui.QMessageBox.question(self, "Deleting minima", 
                                   "Do you want to delete minima %d with energy %g"%(min1._id, min1.energy), 
                                   QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel)
        if(ret == QtGui.QMessageBox.Ok):
            print "deleting minima"
            print "deleting minimum", min1._id, min1.energy
            self.RemoveMinimum(min1)
            self.system.database.removeMinimum(min1)


#this is currently not used.  it may be used later though
#    def LocalConnect(self):
#        self.local_connect = self.system.create_local_connect()
#        
#        min1 = self.ui.listMinima1.selectedItems()[0].minimum
#        min2 = self.ui.listMinima2.selectedItems()[0].minimum
#        res = self.local_connect.connect(min1, min2)
#        ntriplets = len(res.new_transition_states)
#        
#        path = []
#        for i in range(ntriplets):
#            tsret, m1ret, m2ret = res.new_transition_states[i]
#            local_path = []
#            local_path.append(m1ret[0])
#            local_path.append(tsret.coords)
#            local_path.append(m2ret[0])
#            smoothpath = self.system.smooth_path(local_path)
#            path += list(smoothpath)
#        
#        coords = np.array(path)
#        self.nebcoords = coords
#        self.ui.oglPath.setCoords(coords[0,:], 1)
#        self.ui.oglPath.setCoords(None, 2)
#        self.ui.sliderFrame.setRange(0, coords.shape[0]-1)
    
    def doubleEndedConnect(self):
        return self._doubleEndedConnect(reconnect=False)

    def doubleEndedReConnect(self):
        return self._doubleEndedConnect(reconnect=True)

    def _doubleEndedConnect(self, reconnect=False, min1min2=None):
#        min1 = self.ui.listMinima1.selectedItems()[0].minimum
#        min2 = self.ui.listMinima2.selectedItems()[0].minimum
#        min1 = self.ui.
        if min1min2 is None:
            min1 = self.ui.oglPath.minima[1]
            min2 = self.ui.oglPath.minima[2]
        else:
            min1, min2 = min1min2
        database = self.system.database
        double_ended_connect = self.system.get_double_ended_connect(min1, min2, database, 
                                                                       fresh_connect=reconnect)
        double_ended_connect.connect()
        mints, S, energies = double_ended_connect.returnPath()
        clist = [m.coords for m in mints]
        print "done finding path, now just smoothing path.  This can take a while"
        smoothpath = self.system.smooth_path(clist)
        print "done"
        
        coords = np.array(smoothpath)
        self.nebcoords = coords
        self.nebenergies = np.array(energies)
        self.ui.oglPath.setCoords(coords[0,:], 1)
        self.ui.oglPath.setCoords(None, 2)
        self.ui.sliderFrame.setRange(0, coords.shape[0]-1)
        
        if self.usepymol:
            self.pymolviewer.update_coords(self.nebcoords, index=1, delete_all=True)
     
    def connect_in_optim(self):
        """spawn an OPTIM job and retrieve the minima and transition states 
        it finds"""
        min1 = self.ui.oglPath.minima[1]
        min2 = self.ui.oglPath.minima[2]
        existing_minima = set(self.system.database.minima())
        spawner = self.system.get_optim_spawner(min1.coords, min2.coords)
        spawner.run()
        db = self.system.database
        newminima, newts = spawner.load_results(self.system.database)
        for m in newminima:
            if m not in existing_minima:
                self.NewMinimum(m)
            
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

    def merge_minima(self):
        min1 = self.ui.oglPath.minima[1]
        min2 = self.ui.oglPath.minima[2]
        self._merge_minima(min1, min2)
        
        
#def refresh_pl():
    #pl.pause(0.000001)    
    
def run_gui(systemtype, db=None):
    app = QtGui.QApplication(sys.argv)
    
    myapp = MyForm(app, systemtype)
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
#    myapp = MyForm(systemtype)
#    refresh_timer = QtCore.QTimer()
#    refresh_timer.timeout.connect(refresh_pl)
#    refresh_timer.start(0.)
#    
#    myapp.show()
#    sys.exit(app.exec_())
