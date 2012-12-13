from PyQt4 import QtGui
from ljsystem import NewLJDialog

from pygmin.systems import BLJCluster, NotImplemented
 
class BLJSystem(BLJCluster):
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        super(BLJSystem, self).__init__(self.natoms)
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")

    
    def findTS(self, coords):
        raise NotImplemented    
        
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(BLJSystem)
