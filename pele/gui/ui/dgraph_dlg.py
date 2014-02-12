import sys

from PyQt4 import QtGui
from PyQt4.QtGui import QApplication, QWidget, QColorDialog, QInputDialog
from PyQt4.QtCore import pyqtSlot

import matplotlib.colors as col

import dgraph_browser
from pele.utils.disconnectivity_graph import DisconnectivityGraph, database2graph
from pele.storage import Database, TransitionState
from pele.utils.events import Signal
import networkx as nx
import PyQt4


class LabelMinimumAction(QtGui.QAction):
    def __init__(self, minimum, parent=None):
        QtGui.QAction.__init__(self, "add label", parent)
        self.parent = parent
        self.minimum = minimum
        self.triggered.connect(self.__call__)

    def __call__(self, val):
        dialog = QInputDialog(parent=self.parent)
#         dialog.setLabelText("")
        dialog.setLabelText("set label for minimum: " + str(self.minimum.energy))
        dialog.setInputMode(0)
        dialog.exec_()
        if dialog.result():
            label = dialog.textValue()
            self.parent._minima_labels[self.minimum] = label



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
        
        self._tree_cid = None
        self._minima_cid = None
        
        self._minima_labels = dict()

#        # populate the dropdown list with the color names
#        self._colors = sorted(col.cnames.keys())
#        self.ui.comboBox_colour.addItems(self._colors)
#        [self.ui.comboBox_colour.addItem(s) for s in self._colors]
#        self.ui.comboBox_colour.activated[str].connect(self._color_tree)


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
            line = "self.ui.lineEdit_%s.setText(str(%s))" % (keyword, str(v))
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
        
#         self.line_width = 0.5
        self._set_lineEdit("linewidth",  default=0.5)


    def _get_input_parameters(self):
        self.params = self.input_params.copy()
#        self.params.pop("show_minima")
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

        line_width = self.ui.lineEdit_linewidth.text()
        if len(line_width) > 0:
            self.line_width = float(line_width)
        
        self.title = self.ui.lineEdit_title.text()
        

        params["center_gmin"] = self.ui.chkbx_center_gmin.isChecked()
        self.show_minima = self.ui.chkbx_show_minima.isChecked()
        params["order_by_energy"] = self.ui.chkbx_order_by_energy.isChecked()
        params["order_by_basin_size"] = self.ui.chkbx_order_by_basin_size.isChecked()
        params["include_gmin"] = self.ui.chkbx_include_gmin.isChecked()
        self.show_trees = self.ui.chkbx_show_trees.isChecked()

    
#    @pyqtSlot(str)
#    def _color_tree(self, colour):
#        if self.tree_selected is not None: 
#            c = col.hex2color(col.cnames[str(colour)])
#            print "coloring tree", colour, self.tree_selected
#
#            for tree in self.tree_selected.get_all_trees():
#                tree.data["colour"] = c
#            
#            self.redraw_disconnectivity_graph()
##            self.tree_selected = None

    @pyqtSlot()
    def on_btnRedraw_clicked(self):
        self.redraw_disconnectivity_graph()

    @pyqtSlot()
    def on_btnRebuild_clicked(self):
        self.rebuild_disconnectivity_graph()
    

    def redraw_disconnectivity_graph(self):        
        self.params = self._get_input_parameters()
        self._draw_disconnectivity_graph(self.show_minima, self.show_trees)

    def rebuild_disconnectivity_graph(self):        
        self._get_input_parameters()
        self._minima_labels = dict()
        self._build_disconnectivity_graph(**self.params)
        self._draw_disconnectivity_graph(self.show_minima, self.show_trees)

    def _build_disconnectivity_graph(self, **params):
        if self.database is None:
            graph = self.graph
        else:
            db = self.database
            apply_Emax = "Emax" in params and "T" not in params
            if apply_Emax:
                graph = database2graph(db, Emax=params['Emax'])
            else:
                graph = database2graph(db)
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

    def _on_pick_tree(self, event):
        """a matplotlib callback function for when a tree is clicked on"""
        if event.artist != self._treepoints:
#                print "you clicked on something other than a node"
            return True
        ind = event.ind[0]
        self.tree_selected = self._tree_list[ind]
        print "tree clicked on", self.tree_selected
        
        # launch a color selector dialog and color 
        # all subtrees by the selected color
        color_dialog = QColorDialog(parent=self)
        color_dialog.exec_()
        if color_dialog.result():
            color = color_dialog.selectedColor()
            rgba = color.getRgbF() # red green blue alpha
            print "color", rgba
            rgb = rgba[:3]
            for tree in self.tree_selected.get_all_trees():
                tree.data["colour"] = rgb
            
            self.redraw_disconnectivity_graph()
    
    def _on_left_click_minimum(self, minimum):
        print "you clicked on minimum with id", minimum._id, "and energy", minimum.energy
        self.minimum_selected(minimum)
    
    def _on_right_click_minimum(self, minimum):
        menu = QtGui.QMenu("list menu", parent=self)
        
        action1 = LabelMinimumAction(minimum, parent=self)
        menu.addAction(action1)
        
        print "launching menu"
        menu.exec_(QtGui.QCursor.pos())
        
#        dialog = QInputDialog(parent=self)
##         dialog.setLabelText("")
#        dialog.setLabelText("set label for minimum: " + str(minimum.energy))
#        dialog.setInputMode(0)
#        dialog.exec_()
#        if dialog.result():
#            label = dialog.textValue()
#            self._minima_labels[minimum] = label
            
    
    def _on_pick_minimum(self, event):
        """matplotlib event called when a minimum is clicked on"""
        if event.artist != self._minima_points:
#                print "you clicked on something other than a node"
            return True
        ind = event.ind[0]
        min1 = self._minima_list[ind]
        if event.mouseevent.button == 3:
            self._on_right_click_minimum(min1)
        else:
            self._on_left_click_minimum(min1)
 
    def _draw_disconnectivity_graph(self, show_minima=True, show_trees=False):
        ax = self.canvas.axes
        ax.clear()
        ax.hold(True)

        dg = self.dg

        # plot the lines and set up the rest of the plot using the built in function
        # this might change some of the minima x positions, so this has to go before
        # anything dependent on those positions
        dg.plot(axes=ax, show_minima=False, linewidth=self.line_width, 
                title=self.title)
        if len(self._minima_labels) > 0:
            dg.label_minima(self._minima_labels, axes=ax)
            self.ui.widget.canvas.fig.tight_layout()

#        if show_trees
        if self._tree_cid is not None:
            self.canvas.mpl_disconnect(self._tree_cid)
            self._tree_cid = None
        if show_trees:
            # draw the nodes
            tree_list, x_pos, energies = self._get_tree_layout(dg.tree_graph)

            treepoints = ax.scatter(x_pos, energies, picker=5, color='red', alpha=0.5)
            self._treepoints = treepoints
            self._tree_list = tree_list

#            def on_pick_tree(event):
#                if event.artist != treepoints:
#    #                print "you clicked on something other than a node"
#                    return True
#                ind = event.ind[0]
#                self.tree_selected = tree_list[ind]
#                print "tree clicked on", self.tree_selected
#                
#                color_dialog = QColorDialog(parent=self)
#                color_dialog.exec_()
#                color = color_dialog.selectedColor()
#                rgba = color.getRgbF() # red green blue alpha
#                print "color", rgba
#                rgb = rgba[:3]
#                for tree in self.tree_selected.get_all_trees():
#                    tree.data["colour"] = rgb


            self._tree_cid = self.canvas.mpl_connect('pick_event', self._on_pick_tree)


        #draw minima as points and make them interactive
        if self._minima_cid is not None:
            self.canvas.mpl_disconnect(self._minima_cid)
            self._minima_cid = None
        if show_minima:
            xpos, minima = dg.get_minima_layout()
            energies = [m.energy for m in minima]
            self._minima_points = ax.scatter(xpos, energies, picker=5)
            self._minima_list = minima
            
#             def on_pick_min(event):
#                 if event.artist != points:
#     #                print "you clicked on something other than a node"
#                     return True
#                 ind = event.ind[0]
#                 min1 = minima[ind]
#                 print "you clicked on minimum with id", min1._id, "and energy", min1.energy
#                 self.minimum_selected(min1)
            self._minima_cid = self.canvas.mpl_connect('pick_event', self._on_pick_minimum)

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
    from pele.storage.database import Minimum
    g = nx.Graph()
    # js850> It's not strictly necessary to add the minima explicitly here,
    # but for some reason it is much faster if you do (factor of 2).  Even 
    # if this means there are many more minima in the graph.  I'm not sure 
    # why this is.  This step is already often the bottleneck of the d-graph 
    # calculation.
    minima = db.session.query(Minimum).filter(Minimum.energy <= Emax)
    g.add_nodes_from(minima)
    # if we order by energy first and add the transition states with the largest
    # the we will take the smallest energy transition state in the case of duplicates
    ts = db.session.query(TransitionState).filter(TransitionState.energy <= Emax)\
                                          .order_by(-TransitionState.energy)
    for t in ts: 
        g.add_edge(t.minimum1, t.minimum2, ts=t)
    return g


if __name__ == "__main__":
    
    db = Database("lj31.db", createdb=False)
    if len(db.minima()) < 2:
        raise Exception("database has no minima")
    
    app = QApplication(sys.argv)        
    md = DGraphDialog(db)
    md.show()
    md.rebuild_disconnectivity_graph()
    
    sys.exit(app.exec_()) 
        
