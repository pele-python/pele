import sys

from PyQt4 import QtGui
from PyQt4.QtGui import QApplication, QWidget
from PyQt4.QtCore import pyqtSlot

import matplotlib.colors as col

import dgraph_browser
from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
from pygmin.landscape import TSGraph
from pygmin.storage import Database, TransitionState
from pygmin.utils.events import Signal
import networkx as nx


class DGraphWidget(QWidget):
    """
    dialog for showing and modifying the disconnectivity graph
    
    Parameters
    ----------
    database : Database object
    graph : networkx Graph, optional
        you can bypass the database and pass a graph directly.  if you pass the graph,
        pass None as the database
    params : dict
        initialize the values for the disconnectivity graph
    """
    def __init__(self, database, graph=None, params={}, parent=None):
        super(DGraphWidget, self).__init__(parent=parent)
        
        self.database = database
        self.graph = graph
        
        self.ui = dgraph_browser.Ui_Form()
        self.ui.setupUi(self)
        self.canvas = self.ui.widget.canvas
#        self.ui.wgt_mpl_toolbar = NavigationToolbar()
#        self.toolbar = self.
        
        self.input_params = params.copy()
        self.params = {}
        self.set_defaults()
        
        self.minimum_selected = Signal()
        # self.minimum_selected(minim)

#        self.rebuild_disconnectivity_graph()
        self.colour_tree = []
        self.tree_selected = None

        # populate the dropdown list with the color names
        self._colors = sorted(col.cnames.keys())
        self.ui.comboBox_colour.addItems(self._colors)
        [self.ui.comboBox_colour.addItem(s) for s in self._colors]
        self.ui.comboBox_colour.activated[str].connect(self._color_tree)


    def _set_checked(self, keyword, default):
        """utility to set the default values for check boxes
        
        objects must have the name chkbx_keyword
        """
        if keyword in self.input_params:
            v = self.input_params[keyword]
        else:
            v = default
        line = "self.ui.chkbx_%s.setChecked(bool(%d))" % (keyword, v)
        exec(line)

    def _set_lineEdit(self, keyword, default=None):
        """utility to set the default values for lineEdit objects
                
        objects must have the name lineEdit_keyword
        """

        if keyword in self.input_params:
            v = self.input_params[keyword]
        else:
            v = default
        if v is not None:
            line = "self.ui.lineEdit_%s.setText(str(%d))" % (keyword, v)
            exec(line)



    def set_defaults(self):
        self._set_checked("center_gmin", True)
        self._set_checked("show_minima", True)
        self._set_checked("order_by_energy", False)
        self._set_checked("order_by_basin_size", True)
        self._set_checked("include_gmin", True)
        self._set_checked("show_trees", False)
#        self.ui.chkbx_show_minima.setChecked(True)
#        self.ui.chkbx_order_by_energy.setChecked(False)
#        self.ui.chkbx_order_by_basin_size.setChecked(True)
#        self.ui.chkbx_include_gmin.setChecked(True)

        self._set_lineEdit("Emax")
        self._set_lineEdit("subgraph_size")
        self._set_lineEdit("nlevels")


    def _get_input_parameters(self):
        self.params = {}
        params = self.params
        
        Emax = self.ui.lineEdit_Emax.text()
        if len(Emax) > 0:
            self.params["Emax"] = float(Emax)

        subgraph_size = self.ui.lineEdit_subgraph_size.text()
        if len(subgraph_size) > 0:
            self.params["subgraph_size"] = int(subgraph_size)

        nlevels = self.ui.lineEdit_nlevels.text()
        if len(nlevels) > 0:
            self.params["nlevels"] = int(nlevels)

        offset = self.ui.lineEdit_offset.text()
        if len(offset) > 0:
            params["node_offset"] = float(offset)


        params["center_gmin"] = self.ui.chkbx_center_gmin.isChecked()
        self.show_minima = self.ui.chkbx_show_minima.isChecked()
        params["order_by_energy"] = self.ui.chkbx_order_by_energy.isChecked()
        params["order_by_basin_size"] = self.ui.chkbx_order_by_basin_size.isChecked()
        params["include_gmin"] = self.ui.chkbx_include_gmin.isChecked()
        self.show_trees = self.ui.chkbx_show_trees.isChecked()

    
    @pyqtSlot(str)
    def _color_tree(self, colour):
        if self.tree_selected is not None: 
            c = col.hex2color(col.cnames[str(colour)])
            print "coloring tree", colour, self.tree_selected

            for tree in self.tree_selected.get_all_trees():
                tree.data["colour"] = c
            
            self.redraw_disconnectivity_graph()
#            self.tree_selected = None

    @pyqtSlot()
    def on_btnRedraw_clicked(self):
#        if clicked is None: return
        self.redraw_disconnectivity_graph()
        

    def redraw_disconnectivity_graph(self):        
        self.params = self._get_input_parameters()
        self._draw_disconnectivity_graph(self.show_minima, self.show_trees)

    def rebuild_disconnectivity_graph(self):        
        self._get_input_parameters()
        self._build_disconnectivity_graph(**self.params)
        self._draw_disconnectivity_graph(self.show_minima, self.show_trees)

    def _build_disconnectivity_graph(self, **params):
        if self.database is None:
            graph = self.graph
        else:
            db = self.database
            if self.params.has_key('Emax'):
                graph = reduced_db2graph(db, params['Emax'])
            else:
                graphwrapper = TSGraph(db)
                graph = graphwrapper.graph
        dg = DisconnectivityGraph(graph, **params)
        dg.calculate()
        self.dg = dg

    def _get_tree_layout(self, tree):
        treelist = []
        xlist = []
        energies = []
        for tree in tree.get_all_trees():
            xlist.append(tree.data["x"])
            treelist.append(tree)
            if tree.is_leaf():
                energies.append(tree.data["minimum"].energy)
            else:
                energies.append(tree.data["ethresh"])
        return treelist, xlist, energies

    def _draw_disconnectivity_graph(self, show_minima=True, show_trees=False):
        ax = self.canvas.axes
        ax.clear()
        ax.hold(True)

        dg = self.dg

        # plot the lines and set up the rest of the plot using the built in function
        # this might change some of the minima x positions, so this has to go before
        # anything dependent on those positions
        dg.plot(axes=ax, show_minima=False)

#        if show_trees
        if show_trees:
            # draw the nodes
            tree_list, x_pos, energies = self._get_tree_layout(dg.tree_graph)

            treepoints = ax.scatter(x_pos, energies, picker=5, color='red', alpha=0.5)

            def on_pick_tree(event):
                if event.artist != treepoints:
    #                print "you clicked on something other than a node"
                    return True
                ind = event.ind[0]
                self.tree_selected = tree_list[ind]
                print "tree clicked on", self.tree_selected


            self.canvas.mpl_connect('pick_event', on_pick_tree)


        #draw minima as points and make them interactive
        if show_minima:
            xpos, minima = dg.get_minima_layout()
            energies = [m.energy for m in minima]
            points = ax.scatter(xpos, energies, picker=5)        
            
            def on_pick_min(event):
                if event.artist != points:
    #                print "you clicked on something other than a node"
                    return True
                ind = event.ind[0]
                min1 = minima[ind]
                print "you clicked on minimum with id", min1._id, "and energy", min1.energy
                self.minimum_selected(min1)
            self.canvas.mpl_connect('pick_event', on_pick_min)

        self.canvas.draw()


class DGraphDialog(QtGui.QMainWindow):
    def __init__(self, database, graph=None, params={}, parent=None, app=None):
        super(DGraphDialog, self).__init__( parent=parent)
        self.dgraph_widget = DGraphWidget(database, graph, params, parent=self)
        self.setCentralWidget(self.dgraph_widget)
    
    def rebuild_disconnectivity_graph(self):
        self.dgraph_widget.rebuild_disconnectivity_graph()
        

def reduced_db2graph(db, Emax):
    '''
    make a networkx graph from a database including only transition states with energy < Emax
    '''
    g = nx.Graph()
    ts = db.session.query(TransitionState).filter(TransitionState.energy <= Emax).all()
    for t in ts: 
        g.add_edge(t.minimum1, t.minimum2, ts=t)
    return g


if __name__ == "__main__":
    
    db = Database("lj31.db")
    if len(db.minima()) < 2:
        raise Exception("database has no minima")
    
    app = QApplication(sys.argv)        
    md = DGraphDialog(db)
    md.show()
    md.rebuild_disconnectivity_graph()
    
    sys.exit(app.exec_()) 
        
