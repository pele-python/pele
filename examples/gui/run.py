from PyQt4 import QtCore, QtGui
import MainWindow 
import sys
from storage import savenlowest

class QMinimumInList(QtGui.QListWidgetItem):
    def setCoords(self, coords):
        self.coords = coords

class MyForm(QtGui.QMainWindow):
    def __init__(self, storage, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.widget.setCoords(storage.data[0][1])
        self.storage = storage

    def SelectMinimum(self, item):
        self.ui.widget.setCoords(item.coords)
        self.ui.widget.repaint()
        

def run_bh(storage):
    import numpy as np
    import potentials.lj as lj
    import basinhopping as bh
    import take_step.random_displacement as random_displacement

    natoms = 12

    # random initial coordinates
    coords = np.random.random(3 * natoms)

    potential = lj.LJ()#1.0, 1.0, None)

    step = random_displacement.takeStep(stepsize=0.5)

    opt = bh.BasinHopping(coords,potential,
                          temperature=5., takeStep=step.takeStep, storage=storage.insert)
    opt.run(20)

if __name__ == "__main__":
    storage = savenlowest.SaveN(nsave=100)
    run_bh(storage)
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm(storage)
    
    c=1
    for i in storage.data:
        item = QMinimumInList(str(i[0]))
        item.setCoords(i[1].reshape(i[1].size/3, 3))
        myapp.ui.listWidget.insertItem(c, item)
        c+=1
        
    myapp.show()
    sys.exit(app.exec_())
