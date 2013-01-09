import tempfile
import numpy as np

from pygmin.systems import AtomicCluster
from pygmin.potentials import LJ
from pygmin.utils.xyz import write_xyz

__all__ = ["LJCluster"]

class LJCluster(AtomicCluster):
    """
    define the System class for a Lennard-Jones cluster
    
    Parameters
    ----------
    natoms : int 
    
    See Also
    --------
    BaseSystem, AtomicCluster
    """
    def __init__(self, natoms):
        super(LJCluster, self).__init__()
        self.natoms = natoms
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
    
    def get_permlist(self):
        return [range(self.natoms)]
    
    def get_potential(self):
        return LJ()

    #
    #below here is stuff only for the gui
    #

    def draw(self, coordslinear, index):
        from OpenGL import GL,GLUT
        coords = coordslinear.reshape(coordslinear.size/3, 3)
        com=np.mean(coords, axis=0)                  
        for xx in coords:
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.5,30,30)
            GL.glPopMatrix()
    
    def load_coords_pymol(self, coordslist, oname, index=1):
        """load the coords into pymol
        
        the new object must be named self.oname so we can manipulate it later
        
        this will be different for every system
        
        this will create a temporary xyz file and tell pymol to load it from there
        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object will be named oname
        index : int
            which object are we replacing.  can be 1 or 2
        """
        import pymol

        suffix = ".xyz"
    #    f, fname = self.get_tempfile(suffix=suffix)
        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        fname = f.name
        
        
        print "saving files to", fname
        
        from pygmin.mindist import CoMToOrigin
        
        for coords in coordslist:
            coords = CoMToOrigin(coords.copy())
            write_xyz(f, coords, ["LA"])
        f.flush()
        
    #    self.f = f # so f won't be garbage collected and deleted
                
        #load the molecule from the temporary file
        pymol.cmd.load(fname)
        
        #rename the molecule object
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        print "objectname", objectname
        pymol.cmd.set_name(objectname, oname)
        
        
        pymol.cmd.hide("everything", oname)
        pymol.cmd.show("spheres", oname)
        
        if index == 1:
            pymol.cmd.color("red", oname)
        else:
            pymol.cmd.color("gray", oname)


#
#only for testing below here
#

def run():
    #create the system object
    sys = LJCluster(15)
    
    #create a database
    db = sys.create_database()
    
    #do a short basinhopping run
    bh = sys.get_basinhopping(database=db, outstream=None)
    while len(db.minima()) < 2:
        bh.run(100)
    
    #try to connect the lowest two minima
    min1, min2 = db.minima()[:2]
    connect = sys.get_double_ended_connect(min1, min2, db)
    connect.connect()
    

if __name__ == "__main__":
    run()