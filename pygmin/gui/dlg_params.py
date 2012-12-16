import sqlalchemy.orm
from PyQt4 import QtGui
from PyQt4.Qt import QVariant

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
        self.model.itemChanged.connect(self.itemChanged)
        
    def fill(self, node, params):
        i=0
        for key,value in params.iteritems():         
            item = QtGui.QStandardItem(key)
            item.setEditable(False)
            if hasattr(value, "iteritems"):
                self.fill(item, value)
            else:
                #item.setChild(0, 0, QtGui.QStandardItem(str(value)))
                editable = QtGui.QStandardItem(str(value))
                node.setChild(i, 1, editable)
                editable.setData((params, key))
                
            node.setChild(i, 0, item)
            i+=1
    
    def itemChanged(self, item):
        tmp = item.data().toPyObject()
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
    d = {"str": "hello", "subitem": { "int": 1}, "float": 1.0} 
    app = QtGui.QApplication(sys.argv)
    dlg = DlgParams(d)
    dlg.show()
    app.exec_()