import networkx as nx
import numpy as np
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QDialog, QWidget

#from pele.gui.ui.mplwidget import MPLWidgetWithToolbar
from pele.gui.ui.graph_view_ui import Ui_Form
from pele.utils.events import Signal
from pele.utils.disconnectivity_graph import database2graph

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class GraphViewWidget(QWidget):
    def __init__(self, database=None, parent=None, app=None, minima=None):
        QWidget.__init__(self, parent=parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
#        self.widget = GraphDisplayWidget(parent=parent)
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
    
    
    def on_btn_show_all_clicked(self, clicked=None):
        if clicked is None: return
        self.show_all()
        
    def show_all(self):
        self.ui.label_status.setText("showing full graph")
        self.from_minima.clear()
        self.positions.clear()
        self.boundary_nodes.clear()
        self.make_graph()
        self.show_graph()
    
    def make_graph_from(self, minima, cutoff=1):
#        cutoff += 1
        self.from_minima.update(minima)
        minima = self.from_minima
        nodes = set()
        # make a graph from the minima in self.minima and nearest neighbors
        outer_layer = set()
        print "remaking graph"
        for m in minima:
            nodesdir = nx.single_source_shortest_path(self.full_graph, m, cutoff=cutoff)
            print nodesdir
            for n, path in nodesdir.iteritems():
                d = len(path) - 1
                print d, cutoff
                if d < cutoff:
                    print "d <= cutoff", d, cutoff
                    outer_layer.discard(n)
                elif d == cutoff:
                    print "d == cutoff"
                    if n not in nodes:
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
    
    def on_minima_picked_local(self, min1):
        if self.ui.checkBox_zoom.isChecked():
            self.make_graph_from([min1])
            text = "showing graph near minima "
            for m in self.from_minima: 
                text += " " + str(m._id)
            self.ui.label_status.setText(text)
            self.show_graph()
    
    def show_graph(self, fixed=False, show_ids=True):
        import pylab as pl
        if not hasattr(self, "graph"):
            self.make_graph()
        
        print "showing graph"
        ax = self.axes
        ax.clear()
        graph = self.graph
        
        #get the layout of the nodes from networkx
        oldlayout = self.positions
        layout = nx.spring_layout(graph, pos=oldlayout)#, fixed=fixed)
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
        #color the nodes by energy
        e = np.array([m.energy for m, xy in layoutlist])
        #plot the nodes
        points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                            s=markersize, c=e, cmap=pl.cm.autumn)
        if not hasattr(self, "colorbar"):
            self.colorbar = self.fig.colorbar(points)
        
        boundary_points = None
        if self.boundary_nodes:
            # draw the boundary nodes as empty circles with thin lines
            boundary_layout_list = filter(lambda nxy: nxy[0] in self.boundary_nodes, layout.items())
            xypos = np.array([xy for n, xy in boundary_layout_list])
            #plot the nodes
            import matplotlib as mpl
#            marker = mpl.markers.MarkerStyle("o", fillstyle="none")
#            marker.set_fillstyle("none")
            boundary_points = ax.scatter(xypos[:,0], xypos[:,1], picker=5, 
                                s=markersize, marker="o", facecolors="none", linewidths=.5)

        
        #scale the axes so the points are not cutoff
        xmax = max((x for x,y in layout.itervalues() ))
        xmin = min((x for x,y in layout.itervalues() ))
        ymax = max((y for x,y in layout.itervalues() ))
        ymin = min((y for x,y in layout.itervalues() ))
        dx = (xmax - xmin)*.1
        dy = (ymax - ymin)*.1
        ax.set_xlim([xmin-dx, xmax+dx])
        ax.set_ylim([ymin-dy, ymax+dy])
        import matplotlib.pyplot as plt
#        self.fig.relim()
#        self.fig.tight_layout()
#        plt.tight_layout()
#        self.fig.set_tight_layout(True)
        
        
        def on_pick(event):
            ind = event.ind[0]
            if event.artist == points:
                min1 = layoutlist[ind][0]
            elif event.artist == boundary_points:
                min1 = boundary_layout_list[ind][0]
            else:
                return True
            print "you clicked on minimum with id", min1._id, "and energy", min1.energy
            self.on_minima_picked(min1)
            self.on_minima_picked_local(min1)

        
        self.fig.canvas.mpl_connect('pick_event', on_pick)

        
        self.canvas.draw()
        self.app.processEvents()
        

class GraphViewDialog(QtGui.QMainWindow):
    def __init__(self, database, parent=None, app=None):
        QtGui.QMainWindow.__init__(self, parent=parent)

        self.widget = GraphViewWidget(database=database, parent=self, app=app)
        self.setCentralWidget(self.widget)
        
        self.app = app


    def start(self):
        self.widget.show_all()
#        self.widget.make_graph()
#        gmin = self.widget.database.minima()[0:1]
#        self.widget.make_graph_from(gmin)
#        self.widget.show_graph()

def start():
    wnd.start()


if __name__ == "__main__":
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
    
        
    
    
    wnd = GraphViewDialog(db, app=app)
#    decrunner = DECRunner(system, db, min1, min2, outstream=wnd.textEdit_writer)
    glutInit()
    wnd.show()
    from PyQt4.QtCore import QTimer
    QTimer.singleShot(10, start)
    sys.exit(app.exec_()) 
