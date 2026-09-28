from __future__ import print_function
import tempfile
import multiprocessing as mp

from . import pymol


from pele.utils.xyz import write_xyz


class PymolViewer(object):
    """
    this class sets up and maintains a pymol viewer for the gui
    
    Parameters
    ----------
    conn : multiprocessing pipe
        the child end of the commumication pipe
    """
    def __init__(self, load_coords_pymol):
#        mp.Process.__init__(self)
#        self.conn = conn
        
#        import pymol
        print("finishing launching pymol")
        pymol.finish_launching() #  must do this before anything else
        print("done launching pymol")
        
        self.load_coords_pymol = load_coords_pymol
        
        self.oname1 = "molecule1" #  name of the object
        self.oname2 = "molecule2" #  name of the object
        
    
#    def get_tempfile(self, suffix=".xyz"):
#        """create a temporary file and return that open file and the file name
#        
#        the file will be delete when f is garbage collected
#        """
#        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
#        return f, f.name
            
    def update_coords(self, coordslist, index=1, delete_all=False):
        if index == 1: oname = self.oname1
        else : oname = self.oname2

        if delete_all:
            pymol.cmd.delete("all")
        else:
            pymol.cmd.delete(oname)
        self.load_coords_pymol(coordslist, oname, index=index)
        
    def __del__(self):
        print("quitting pymol safely")
        pymol.cmd.quit()
        
        

        
