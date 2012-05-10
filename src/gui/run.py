from PyQt4 import QtCore, QtGui, Qt
import PyQt4
import MainWindow 
import sys
import bhrunner

class QMinimumInList(QtGui.QListWidgetItem):
    def setCoords(self, coords):
        self.coords = coords
    def setMinimumId(self, id):
        self.minid = id

class MyForm(QtGui.QMainWindow):
    def __init__(self, systemtype, parent=None):
        QtGui.QWidget.__init__(self)
        self.ui = MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
	self.systemtype = systemtype
        self.NewSystem()
        
    def NewSystem(self):
        self.system = self.systemtype()
        self.system.storage.onMinimumAdded=self.NewMinimum
        self.system.storage.onMinimumRemoved=self.RemoveMinimum
    
    def SelectMinimum(self, item):
        self.ui.widget.setSystem(self.system)
        self.ui.widget.setCoords(item.coords)
        self.ui.widget.repaint()
        
    def NewMinimum(self, minimum):
        E=minimum[0]
        id=minimum[2]
        coords=minimum[1]
        c=1
        #myapp.ui.listWidget.clear()
        #for i in self.bhrunner.storage.data:
        #    E = i[0]
        #    coords = i[1]
        item = QMinimumInList('%.4f'%E)
        item.setCoords(coords.reshape(coords.size/3, 3))
        item.setMinimumId(id)
        #myapp.ui.listWidget.insertItem(1, item)
        self.ui.listWidget.addItem(item)
        #myapp.ui.listWidget.insertItem(1, str(E))
        #    c=c+1
        self.ui.listWidget.sortItems(1) #Qt.AscendingOrder)
    
    def RemoveMinimum(self, min):
        E = min[0]
        id = min[2]
        widget = self.ui.listWidget
        itms = widget.findItems('*', QtCore.Qt.MatchWildcard)
        for i in itms:
            if(i.minid == id):
                widget.takeItem(widget.row(i)) 
        
    def StartBasinHopping(self):
        self.bhrunner = bhrunner.BHRunner(self.system)
        self.bhrunner.start()
    
def run_gui(systemtype):
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm(systemtype)
            
    myapp.show()
    sys.exit(app.exec_())
