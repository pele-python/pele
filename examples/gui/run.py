from PyQt4 import QtCore, QtGui, Qt
import PyQt4
import MainWindow 
import sys
from storage import savenlowest
import bhrunner

class QMinimumInList(QtGui.QListWidgetItem):
    def setCoords(self, coords):
        self.coords = coords
    def setMinimumId(self, id):
        self.minid = id

class MyForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self)
        self.ui = MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        import ljsystem
        self.system=ljsystem.LJSystem()

    def NewSystem(self):
        import ljsystem
        self.system = ljsystem.LJSystem()
    
    def SelectMinimum(self, item):
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
        myapp.ui.listWidget.addItem(item)
        #myapp.ui.listWidget.insertItem(1, str(E))
        #    c=c+1
        myapp.ui.listWidget.sortItems(1) #Qt.AscendingOrder)
    
    def RemoveMinimum(self, min):
        E = min[0]
        id = min[2]
        widget = myapp.ui.listWidget
        itms = widget.findItems('*', QtCore.Qt.MatchWildcard)
        for i in itms:
            if(i.minid == id):
                widget.takeItem(widget.row(i)) 
        
    def StartBasinHopping(self):
        self.bhrunner = bhrunner.BHRunner(self.system,
                                  onMinimumAdded=self.NewMinimum,
                                  onMinimumRemoved=self.RemoveMinimum)
        self.bhrunner.start()
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm()
            
    myapp.show()
    sys.exit(app.exec_())
