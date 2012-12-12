from PyQt4 import QtGui
import NewLJ
from ljsystem import NewLJDialog

from pygmin.transition_states import NEB, InterpolatedPath

from pygmin.systems import BLJCluster, NotImplemented
 
class BLJSystem(BLJCluster):
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")

        super(BLJSystem, self).__init__(self.natoms)
    
    def findTS(self, coords):
        raise NotImplemented    
        
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(BLJSystem)
