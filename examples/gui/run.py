from PyQt4 import QtCore, QtGui
import MainWindow 
import sys
from storage import savenlowest
import bhrunner

class QMinimumInList(QtGui.QListWidgetItem):
    def setCoords(self, coords):
        self.coords = coords

class MyForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self)
        self.ui = MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.bhrunner = bhrunner.BHRunner(onNewMinimumFound=self.NewMinimum)

    def SelectMinimum(self, item):
        self.ui.widget.setCoords(item.coords)
        self.ui.widget.repaint()
        
    def NewMinimum(self, E, coords):
        print "New minimum",E
        c=1
        #myapp.ui.listWidget.clear()
        #for i in self.bhrunner.storage.data:
        #    E = i[0]
        #    coords = i[1]
        item = QMinimumInList(str(100-E))
        item.setCoords(coords.reshape(coords.size/3, 3))
        #myapp.ui.listWidget.insertItem(1, item)
        myapp.ui.listWidget.addItem(item)
        #myapp.ui.listWidget.insertItem(1, str(E))
        #    c=c+1
        
    def StartBasinHopping(self):
        self.bhrunner.start()
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm()
            
    myapp.show()
    sys.exit(app.exec_())
