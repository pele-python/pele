import sqlalchemy.orm
from PyQt4 import QtGui

from ui_params import Ui_Dialog as UI

class DlgParams(QtGui.QDialog):
    def __init__(self, params):
        QtGui.QDialog.__init__(self)
        self.ui = UI()
        self.ui.setupUi(self)
        self.params = params

        self.model = QtGui.QStandardItemModel()
        self.ui.treeParams.setModel(self.model)
        
        root = QtGui.QStandardItem("pygmin")
        self.model.appendRow(root)
        self.model.setColumnCount(2)
        self.fill(root, params)
    
    def fill(self, node, params):
        i=0
        for key,value in params.iteritems():         
            item = QtGui.QStandardItem(key)
            item.setEditable(False)
            if hasattr(value, "iteritems"):
                self.fill(item, value)
            else:
                #item.setChild(0, 0, QtGui.QStandardItem(str(value)))
                node.setChild(i, 1, QtGui.QStandardItem(str(value)))
                print str(value)
                
            node.setChild(i, 0, item)
            i+=1
            
                
    def accept(self, *args, **kwargs):
        return QtGui.QDialog.accept(self, *args, **kwargs)
    
    def reject(self, *args, **kwargs):
        return QtGui.QDialog.reject(self, *args, **kwargs)


