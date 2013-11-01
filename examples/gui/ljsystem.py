from PyQt4 import QtGui
import NewLJ
import numpy as np

from pele.config import config

#from pele.transition_states import dimer, tstools

from pele.systems import LJCluster
 
class LJSystem(LJCluster):
    def __init__(self):
        dlg = NewLJDialog()
        dlg.exec_()
        self.natoms = dlg.natoms()
        if dlg.result() == QtGui.QDialog.Rejected:
            raise BaseException("Aborted parameter dialog")
        super(LJSystem, self).__init__(self.natoms)
#        self.params.gui.basinhopping_nsteps = 1000
    
#    def get_optim_spawner(self, coords1, coords2):
#        from pele.systems.spawn_OPTIM import SpawnOPTIM_LJ
#        import os
##        # TODO: this should be passable somehow
##        optim = os.path.expanduser("~")+"/git/OPTIM/source/build/OPTIM"
###        optim = "OPTIM"
#        optim = config.get("exec", "OPTIM")
#        optim = os.path.expandvars(os.path.expanduser(optim))
#        print "optim executable", optim
#        return SpawnOPTIM_LJ(coords1, coords2, self, OPTIM=optim, tempdir=True)


class NewLJDialog(QtGui.QDialog,NewLJ.Ui_DialogLJSetup):
    def __init__(self):
        QtGui.QDialog.__init__(self)
        self.setupUi(self)
    def natoms(self):
        return int(self.lineNatoms.text())
    def nsave(self):
        return int(self.lineNsave.text())
        
if __name__ == "__main__":
    import pele.gui.run as gr
    gr.run_gui(LJSystem)
