import sqlalchemy.orm
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import QVariant

from ui_params import Ui_Dialog as UI

class DlgParams(QtGui.QDialog):
    def __init__(self, params):
        QtGui.QDialog.__init__(self)
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
        
    def fill(self, params, node=None):
        i=0
        for key,value in params.iteritems():         
            new_node = QtGui.QStandardItem(str(key))
            new_node.setEditable(False)
            new_node.setData((params, key))
            if hasattr(value, "iteritems"):
                self.fill(value, new_node)
                editable= QtGui.QStandardItem()
                editable.setEditable(False)
                editable.setEnabled(False)
            else:
                #item.setChild(0, 0, QtGui.QStandardItem(str(value)))
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
        if(hasattr(params[key], "iteritems")):
            menu.addAction(self.tr("Add option"))
        else:
            menu.addAction(self.tr("Delete"))

        menu.exec_(self.ui.treeParams.viewport().mapToGlobal(position))
    
    def item_changed(self, item):
        tmp = item.data().toPyObject()
        print tmp
        if tmp is None: return
        dict_, attr_ = tmp
        try:
            dict_[attr_] = type(dict_[attr_])(item.text())
        except ValueError:
            item.setText(str(dict_[attr_]))
             
    def accept(self, *args, **kwargs):
        return QtGui.QDialog.accept(self, *args, **kwargs)
    
    def reject(self, *args, **kwargs):
        return QtGui.QDialog.reject(self, *args, **kwargs)

if __name__ == "__main__":
    import sys
    d = {"str": "hello", "subitem": { "int": 1, "subsub": {"test": 3}}, "float": 1.0} 
    app = QtGui.QApplication(sys.argv)
    dlg = DlgParams(d)
    dlg.show()
    app.exec_()