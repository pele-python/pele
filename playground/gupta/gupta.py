import gmin_ as GMIN
from pele.potentials import GMINPotential
from pele.systems import LJCluster
import numpy as np

class GUPTASystem(LJCluster):
    def __init__(self):
        GMIN.initialize()
        self.natoms = GMIN.getNAtoms()
        super(GUPTASystem, self).__init__(self.natoms)

        qp = self.params.structural_quench_params
        qp["tol"]=1e-5
        qp["maxErise"]=1e-5
        qp["maxstep"]=0.1
        qp["iprint"]=-1
        qp["debug"]=False
        
        neb = self.params.double_ended_connect.local_connect_params.NEBparams
        neb["image_density"]=5
        neb["adjustk_freq"]=5

    
    def get_potential(self):
        return GMINPotential(GMIN)
    
if __name__ == "__main__":
    import pele.gui.run as gr
    gr.run_gui(GUPTASystem)