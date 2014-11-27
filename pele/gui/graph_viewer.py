import networkx as nx
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QWidget

from pele.gui.ui.graph_view_ui import Ui_Form
from pele.utils.events import Signal
from pele.utils.disconnectivity_graph import database2graph
from pele.gui.ui.dgraph_dlg import minimum_energy_path, check_thermodynamic_info
from pele.rates import compute_committors



try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class ShowPathAction(QtGui.QAction):
    """this action will show the minimum energy path to minimum1"""
    def __init__(self, minimum1, minimum2, parent=None):
        QtGui.QAction.__init__(self, "show path to %d" % minimum2._id, parent)
        self.parent = parent
        self.minimum1 = minimum1
        self.minimum2 = minimum2
        self.triggered.connect(self.__call__)

    def __call__(self, val):
        self.parent._show_minimum_energy_path(self.minimum1, self.minimum2)

class ColorByCommittorAction(QtGui.QAction):
    """this action will color the graph by committor probabilities"""
    def __init__(self, minimum1, minimum2, parent=None):
        QtGui.QAction.__init__(self, "color by committor %d" % minimum2._id, parent)
        self.parent = parent
        self.minimum1 = minimum1
        self.minimum2 = minimum2
        self.triggered.connect(self.__call__)

    def __call__(self, val):
        dialog = QtGui.QInputDialog(parent=self.parent)
        dialog.setLabelText("Temperature for committor calculation")
        dialog.setInputMode(2)
        dialog.setDoubleValue(1.)
        dialog.exec_()
        if dialog.result():
            T = dialog.doubleValue()
            self.parent._color_by_committor(self.minimum1, self.minimum2, T=T)



class GraphViewWidget(QWidget):
    def __init__(self, database=None, parent=None, app=None, minima=None):
        QWidget.__init__(self, parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        self.database = database
        self.minima = minima
        self.app = app
        
        self.canvas = self.ui.canvas.canvas
        
        self.axes = self.canvas.axes
        self.fig = self.canvas.fig
        
        self.on_minima_picked = Signal()
        
        self.from_minima = set()
        self.positions = dict()
        self.boundary_nodes = set()
        
        self._selected_minimum = None
        self._mpl_cid = None
        
        self._minima_color_value = None
    
    
    def on_btn_show_all_clicked(self, clicked=None):
        if clicked is None: return
        self.show_all()
    
    def _reset_minima_lists(self):
        self.from_minima.clear()
        self.positions.clear()
        self.boundary_nodes.clear()
        self._minima_color_value = None
    
    def show_all(self):
        self.ui.label_status.setText("showing full graph")
        self._reset_minima_lists()
        self.make_graph()
        self.show_graph()
    
    def make_graph_from(self, minima, cutoff=1):
        """rebuild the graph using only the passed minima and those in self.from_minima"""
#        cutoff += 1
        self.from_minima.update(minima)
        minima = self.from_minima
        nodes = set()
        # make a graph from the minima in self.minima and nearest neighbors
        outer_layer = set()
        for m in minima:
            nodesdir = nx.single_source_shortest_path(self.full_graph, m, cutoff=cutoff)
            for n, path in nodesdir.iteritems():
                d = len(path) - 1
                if d < cutoff:
                    # n is close to m, remove it from outer layer
                    outer_layer.discard(n)
                elif d == cutoff:
                    if n not in nodes:
                        # n is in the outer layer of m and not near any other nodes.
                        outer_layer.add(n)
            nodes.update(nodesdir)

        self.boundary_nodes = outer_layer
        self.graph = self.full_graph.subgraph(nodes)

        # remove nodes not in the graph from the dictionary positions
        difference = set(self.positions.viewkeys())
        difference.difference_update(self.graph.nodes())
        for m in difference:
            self.positions.pop(m)

        print "boundary nodes", len(self.boundary_nodes), self.graph.number_of_nodes()
    
    def make_graph(self, database=None, minima=None):
        """build an nx graph from the database"""
        if database is None:
            database = self.database
        if minima is None:
            minima = self.minima

        print "making graph", database, minima
        # get the graph object, eliminate nodes without edges
        graph = database2graph(database)
        if minima is not None:
            to_remove = set(graph.nodes()).difference(set(minima))
            graph.remove_nodes_from(to_remove)
        self.full_graph = graph
        print graph.number_of_nodes()
        degree = graph.degree()
        nodes = [n for n, nedges in degree.items() if nedges > 0]
        self.graph = graph.subgraph(nodes)
        print self.graph.number_of_nodes(), self.graph.number_of_edges()
    
    def _show_minimum_energy_path(self, m1, m2):
        """show only the minima in the path from m1 to m2"""
        self._reset_minima_lists()
        path = minimum_energy_path(self.full_graph, m1, m2)
        self.make_graph_from(path)
        print "there are", len(path), "minima in the path from", m1._id, "to", m2._id 
        status = "showing path from minimum %d to %d" % (m1._id, m2._id)
        self.ui.label_status.setText(status)
        self.show_graph()
    
    
    def _color_by_committor(self, min1, min2, T=1.):
        print "coloring by the probability that a trajectory gets to minimum", min1._id, "before", min2._id
        # get a list of transition states in the same cluster as min1
        edges = nx.bfs_edges(self.graph, min1)
        transition_states = [ self.graph.get_edge_data(u, v)["ts"] for u, v in edges ]
        if not check_thermodynamic_info(transition_states):
            raise Exception("The thermodynamic information is not yet computed.  Use the main menu of the gui under 'Actions'")
        
        A = [min2]
        B = [min1]
        committors = compute_committors(transition_states, A, B, T=T)
        self._minima_color_value = committors
        def get_committor(m):
            try:
                return committors[m]
            except KeyError:
                return 1.
        self._minima_color_value = get_committor
        self.show_graph()

    
    def _on_right_click_minimum(self, minimum):
        """create a menu with the list of available actions"""
        print "you right clicked on minimum with id", minimum._id, "and energy", minimum.energy
        menu = QtGui.QMenu("list menu", parent=self)
        
        if self._selected_minimum is not None:
            menu.addAction(ShowPathAction(minimum, self._selected_minimum, parent=self))
            menu.addAction(ColorByCommittorAction(minimum, self._selected_minimum, parent=self))

        menu.exec_(QtGui.QCursor.pos())

    def _on_left_click_minimum(self, min1):
        self._selected_minimum = min1
        print "you clicked on minimum with id", min1._id, "and energy", min1.energy
        self.on_minima_picked(min1)
        if self.ui.checkBox_zoom.isChecked():
            self.make_graph_from([min1])
            text = "showing graph near minima "
            for m in self.from_minima: 
                text += " " + str(m._id)
            self.ui.label_status.setText(text)
            self.show_graph()
    
    def _on_mpl_pick_event(self, event):
        """matplotlib event called when a minimum is clicked on"""
        artists = {self._boundary_points, self._minima_points}
        if event.artist not in artists:
#                print "you clicked on something other than a node"
            return True
        ind = event.ind[0]
        if event.artist == self._minima_points:
            min1 = self._mimima_layout_list[ind][0]
        else:
            min1 = self._boundary_layout_list[ind][0]
        if event.mouseevent.button == 3:
            self._on_right_click_minimum(min1)
        else:
            self._on_left_click_minimum(min1)


    
    def show_graph(self, fixed=False, show_ids=True):
        """draw the graph"""
        import pylab as pl
        if not hasattr(self, "graph"):
            self.make_graph()
        
        print "showing graph"
        # I need to clear the figure and make a new axes object because
        # I can't find any other way to remove old colorbars
        self.fig.clf()
        self.axes = self.fig.add_subplot(111)
        ax = self.axes
        ax.clear()
        graph = self.graph
        
        # get the layout of the nodes from networkx
        oldlayout = self.positions
        layout = nx.spring_layout(graph, pos=oldlayout)
        self.positions.update(layout)
        layout = self.positions
        
        # draw the edges as lines
        from matplotlib.collections import LineCollection
        linecollection = LineCollection([(layout[u], layout[v]) for u, v in graph.edges()
                 if u not in self.boundary_nodes and v not in self.boundary_nodes])
        linecollection.set_color('k')
        ax.add_collection(linecollection)

        if self.boundary_nodes:
            # draw the edges connecting the boundary nodes as thin lines
            from matplotlib.collections import LineCollection
            linecollection = LineCollection([(layout[u], layout[v]) for u, v in graph.edges()
                     if u in self.boundary_nodes or v in self.boundary_nodes])
            linecollection.set_color('k')
            linecollection.set_linewidth(0.2)
            ax.add_collection(linecollection)

        markersize = 8**2
        
        # draw the interior nodes
        interior_nodes = set(graph.nodes()) - self.boundary_nodes
        layoutlist = filter(lambda nxy: nxy[0] in interior_nodes, layout.items())
        xypos = np.array([xy for n, xy in layoutlist])
        if self._minima_color_value is None:
            #color the nodes by energy
            color_values = [m.energy for m, xy in layoutlist]
        else:
            color_values = [self._minima_color_value(m) for m, xy in layoutlist]
        #plot the nodes
        self._minima_points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                            s=markersize, c=color_values, cmap=pl.cm.autumn)
        self._mimima_layout_list = layoutlist
#        if not hasattr(self, "colorbar"):
        self.colorbar = self.fig.colorbar(self._minima_points)
        
        self._boundary_points = None
        self._boundary_list = []
        if self.boundary_nodes:
            # draw the boundary nodes as empty circles with thin lines
            boundary_layout_list = filter(lambda nxy: nxy[0] in self.boundary_nodes, layout.items())
            xypos = np.array([xy for n, xy in boundary_layout_list])
            #plot the nodes
#            marker = mpl.markers.MarkerStyle("o", fillstyle="none")
#            marker.set_fillstyle("none")
            self._boundary_points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                                s=markersize, marker="o", facecolors="none", linewidths=.5)
            self._boundary_layout_list = boundary_layout_list

        
        #scale the axes so the points are not cutoff
        xmax = max((x for x,y in layout.itervalues() ))
        xmin = min((x for x,y in layout.itervalues() ))
        ymax = max((y for x,y in layout.itervalues() ))
        ymin = min((y for x,y in layout.itervalues() ))
        dx = (xmax - xmin)*.1
        dy = (ymax - ymin)*.1
        ax.set_xlim([xmin-dx, xmax+dx])
        ax.set_ylim([ymin-dy, ymax+dy])
        
        if self._mpl_cid is not None:
            self.canvas.mpl_disconnect(self._mpl_cid)
            self._mpl_cid = None
        self._mpl_cid = self.fig.canvas.mpl_connect('pick_event', self._on_mpl_pick_event)

        
        self.canvas.draw()
        self.app.processEvents()
        

class GraphViewDialog(QtGui.QMainWindow):
    def __init__(self, database, parent=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)
        self.setWindowTitle("Connectivity graph")

        self.widget = GraphViewWidget(database=database, parent=self, app=app)
        self.setCentralWidget(self.widget)
        
        self.app = app


    def start(self):
        self.widget.show_all()


def test():
    from OpenGL.GLUT import glutInit
    import sys
    import pylab as pl

    app = QtGui.QApplication(sys.argv)
    from pele.systems import LJCluster
    pl.ion()
    natoms = 13
    system = LJCluster(natoms)
    system.params.double_ended_connect.local_connect_params.NEBparams.iter_density = 5.
    dbname = "lj%dtest.db" % (natoms,)
    db = system.create_database(dbname)
    
    #get some minima
    if False:
        bh = system.get_basinhopping(database=db)
        bh.run(10)
        minima = db.minima()
    else:
        x1, e1 = system.get_random_minimized_configuration()[:2]
        x2, e2 = system.get_random_minimized_configuration()[:2]
        min1 = db.addMinimum(e1, x1)
        min2 = db.addMinimum(e2, x2)
        minima = [min1, min2]

    
    # connect some of the minima
    nmax = min(3, len(minima))
    m1 = minima[0]
    for m2 in minima[1:nmax]:
        connect = system.get_double_ended_connect(m1, m2, db)
        connect.connect()
    
        
    if True:
        from pele.thermodynamics import get_thermodynamic_information
        get_thermodynamic_information(system, db, nproc=4)
    
    wnd = GraphViewDialog(db, app=app)
#    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    glutInit()
    wnd.show()
    from PyQt4.QtCore import QTimer
    def start():
        wnd.start()

    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 

if __name__ == "__main__":
    test()
