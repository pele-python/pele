from PyQt4 import QtCore, QtGui
import MainWindow 
import sys
import bhrunner
import copy

class QMinimumInList(QtGui.QListWidgetItem):
    def setCoords(self, coords):
        self.coords = coords
    def setMinimumId(self, minid):
        self.minid = minid

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
        
    def save(self):
        import pickle
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '.')
        output = open(filename, "w")
        pickle.dump(self.system.storage, output)
        self.system.storage.onMinimumAdded=self.NewMinimum
        self.system.storage.onMinimumRemoved=self.RemoveMinimum
       
    def load(self):
        import pickle
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File', '.')
        infile = open(filename, "r")
        self.system.storage = pickle.load(infile)
        for minimum in self.system.storage.data:
            self.NewMinimum(minimum)
            
    def SelectMinimum(self, item):
        self.ui.widget.setSystem(self.system)
        self.ui.widget.setCoords(item.coords)

    def SelectMinimum1(self, item):
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(item.coords, index=1)
        self.neb = None
    
    def SelectMinimum2(self, item):
        self.ui.oglPath.setSystem(self.system)
        self.ui.oglPath.setCoords(item.coords, index=2)
        self.neb = None
        
    def AlignMinima(self):
        coords1 = self.ui.oglPath.coords[1]
        coords2 = self.ui.oglPath.coords[2]
        coords1, coords2 = self.system.Align(coords1, coords2)
        self.ui.oglPath.setCoords(coords1, 1)
        self.ui.oglPath.setCoords(coords2, 2)
        pass    
    
    def ConnectMinima(self):
        self.neb = self.system.createNEB(self.ui.oglPath.coords[1], self.ui.oglPath.coords[2])
        self.neb.optimize()
        self.ui.oglPath.setCoords(self.neb.coords[0,:], 1)
        self.ui.oglPath.setCoords(None, 2)
        self.ui.sliderFrame.setRange(0,self.neb.coords.shape[0]-1)

    def showFrame(self, i):
        if(self.neb):
            self.ui.oglPath.setCoords(self.neb.coords[i,:])
    
    def showEnergies(self):
        import pylab as pl
        neb = self.neb
        pl.plot(neb.energies, "o-", label="neb energies")
        cl=[]
        en=[]
        for i in xrange(len(neb.energies)):
            if(neb.isclimbing[i]):
                print "climbing image :", i, neb.energies[i]
                cl.append(i)
                en.append(neb.energies[i])
                
        pl.plot(cl, en, "s", label="climbing images", markersize=10, markerfacecolor="none", markeredgewidth=2)
        pl.legend(loc='best')
        pl.show()
     
    def NewMinimum(self, minimum):
        self.AddMinimumToList(self.ui.listWidget, minimum)
        self.AddMinimumToList(self.ui.listMinima1, minimum)
        self.AddMinimumToList(self.ui.listMinima2, minimum)

    def AddMinimumToList(self, obj, minimum):
        E=minimum.E
        minid=id(minimum)
        coords=minimum.coords
        item = QMinimumInList('%.4f'%E)
        item.setCoords(coords)
        item.setMinimumId(minid)
        obj.addItem(item)    
        obj.sortItems(1)
        
    
    def RemoveMinimum(self, minimum):
        self.RemoveMinimumFromList(self.ui.listWidget,  minimum)
        self.RemoveMinimumFromList(self.ui.listMinima1,  minimum)
        self.RemoveMinimumFromList(self.ui.listMinima2,  minimum)

    def RemoveMinimumFromList(self, obj, minimum):
        minid = id(minimum)
        itms = obj.findItems('*', QtCore.Qt.MatchWildcard)
        for i in itms:
            if(i.minid == minid):
                obj.takeItem(obj.row(i))
                
    def StartBasinHopping(self):
        self.bhrunner = bhrunner.BHRunner(self.system)
        self.bhrunner.start()
    
def run_gui(systemtype):
    app = QtGui.QApplication(sys.argv)
    myapp = MyForm(systemtype)
            
    myapp.show()
    sys.exit(app.exec_())
