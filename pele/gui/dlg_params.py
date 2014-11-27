from PyQt4 import QtGui, QtCore

from pele.gui.ui_params import Ui_Dialog as UI

class EditParamsWidget(QtGui.QWidget):
    def __init__(self, parent=None, params=None):
        if params is None: params = dict()
        QtGui.QWidget.__init__(self, parent)
        self.params = params
        self.view = QtGui.QTreeView(self)
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.view)
        self.setLayout(vbox)
        
        self.model = QtGui.QStandardItemModel()
        self.model.setHorizontalHeaderLabels(["parameter", "value"])
        self.view.setModel(self.model)
        self.view.setColumnWidth(0, 300)
        self.model.setColumnCount(2)
        
        self.fill(params)
        self.model.itemChanged.connect(self.item_changed)
        self.view.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        
        self.view.setSortingEnabled(True)
        self.view.sortByColumn(0, QtCore.Qt.AscendingOrder)
        
    def fill(self, params, node=None):
        i=0
        for key,value in params.iteritems():
            if callable(value):
                continue
            new_node = QtGui.QStandardItem(str(key))
            new_node.setEditable(False)
            new_node.setData((params, key))
            if hasattr(value, "iteritems"):
                self.fill(value, new_node)
                editable= QtGui.QStandardItem()
                editable.setEditable(False)
                editable.setEnabled(False)
            else:
                if type(value) == bool:
                    editable = QtGui.QStandardItem()
                    editable.setCheckable(True)
                    editable.setEditable(False)
                    editable.setCheckState(QtCore.Qt.Checked if value else QtCore.Qt.Unchecked)
                else:
                    editable = QtGui.QStandardItem(str(value))
                editable.setData((params, key))
                
            if node is None:
                self.model.appendRow([new_node, editable])
            else:
                node.setChild(i, 0, new_node)
                node.setChild(i, 1, editable)
                
            i+=1
            
    def open_context_menu(self, position):
        indexes = self.ui.treeParams.selectedIndexes()
        if len(indexes) == 0:
            return
        d = indexes[0].data(role=QtCore.Qt.UserRole+1).toPyObject()
        if d is None:
            return
            
        params, key = d
        menu = QtGui.QMenu()
        if hasattr(params[key], "iteritems"):
            menu.addAction(self.tr("Add option"))
        else:
            menu.addAction(self.tr("Delete"))

        menu.exec_(self.ui.treeParams.viewport().mapToGlobal(position))
    
    def item_changed(self, item):
        tmp = item.data().toPyObject()
        if tmp is None: return
        dict_, attr_ = tmp
        try:
            if type(dict_[attr_]) == bool:
                dict_[attr_] = True if item.checkState() == QtCore.Qt.Checked else False
            else:
                dict_[attr_] = type(dict_[attr_])(item.text())
        except ValueError:
            item.setText(str(dict_[attr_]))
        
class DlgParams(QtGui.QDialog):
    def __init__(self, params, parent=None):
        QtGui.QDialog.__init__(self, parent=parent)
        self.ui = UI()
        self.ui.setupUi(self)
        self.params = params

        self.model = QtGui.QStandardItemModel()
        self.model.setHorizontalHeaderLabels(["parameter", "value", "description"])
        self.ui.treeParams.setModel(self.model)
        self.ui.treeParams.setColumnWidth(0, 300)
        self.model.setColumnCount(2)
        self.fill(params)
        self.model.itemChanged.connect(self.item_changed)
        self.ui.treeParams.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.treeParams.customContextMenuRequested.connect(self.open_context_menu)
        
        self.ui.treeParams.setSortingEnabled(True)
        self.ui.treeParams.sortByColumn(0, QtCore.Qt.AscendingOrder)
        
    def fill(self, params, node=None):
        i=0
        for key,value in params.iteritems():
            if callable(value):
                continue
                        
            new_node = QtGui.QStandardItem(str(key))
            new_node.setEditable(False)
            new_node.setData((params, key))
            if hasattr(value, "iteritems"):
                self.fill(value, new_node)
                editable= QtGui.QStandardItem()
                editable.setEditable(False)
                editable.setEnabled(False)
            else:
                if type(value) == bool:
                    editable = QtGui.QStandardItem()
                    editable.setCheckable(True)
                    editable.setEditable(False)
                    editable.setCheckState(QtCore.Qt.Checked if value else QtCore.Qt.Unchecked)
                else:
                    editable = QtGui.QStandardItem(str(value))
                editable.setData((params, key))
                
            if node is None:
                self.model.appendRow([new_node, editable])
            else:
                node.setChild(i, 0, new_node)
                node.setChild(i, 1, editable)
                
            i+=1
            
    def open_context_menu(self, position):
        indexes = self.ui.treeParams.selectedIndexes()
        if len(indexes) == 0:
            return
        d = indexes[0].data(role=QtCore.Qt.UserRole+1).toPyObject()
        if d is None:
            return
            
        params, key = d
        menu = QtGui.QMenu()
        if hasattr(params[key], "iteritems"):
            menu.addAction(self.tr("Add option"))
        else:
            menu.addAction(self.tr("Delete"))

        menu.exec_(self.ui.treeParams.viewport().mapToGlobal(position))
    
    def item_changed(self, item):
        tmp = item.data().toPyObject()
        if tmp is None: return
        dict_, attr_ = tmp
        try:
            if type(dict_[attr_]) == bool:
                dict_[attr_] = True if item.checkState() == QtCore.Qt.Checked else False
            else:
                dict_[attr_] = type(dict_[attr_])(item.text())
        except ValueError:
            item.setText(str(dict_[attr_]))
             
    def accept(self, *args, **kwargs):
        return QtGui.QDialog.accept(self, *args, **kwargs)
    
    def reject(self, *args, **kwargs):
        return QtGui.QDialog.reject(self, *args, **kwargs)

if __name__ == "__main__":
    import sys
    d = {"str": "hello", "subitem": { "int": 1, "subsub": {"test": 3}}, "float": 1.0, "bool": True} 
    app = QtGui.QApplication(sys.argv)
    dlg = DlgParams(d)
    dlg.show()
    app.exec_()