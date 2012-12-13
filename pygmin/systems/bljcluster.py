import numpy as np

from pygmin.systems import AtomicCluster
from pygmin.potentials.ljpshiftfast import LJpshift

__all__ = ["BLJCluster"]

class BLJCluster(AtomicCluster):
    """
    define the System class for a lennard jones cluster
    """
    def __init__(self, natoms, ntypeA="default", **potential_kwargs):
        super(BLJCluster, self).__init__()
        self.natoms = natoms
        if ntypeA == "default":
            self.ntypeA = int(self.natoms * 0.8)
        else:
            self.ntypeA = ntypeA
        self.potential_kwargs = potential_kwargs
        

        self.params["database"]["accuracy"] = 1e-3

    
    def get_potential(self):
        return LJpshift(self.natoms, self.ntypeA, **self.potential_kwargs)
    
    def get_permlist(self):
        return [range(self.ntypeA), range(self.ntypeA, self.natoms)]

    def draw(self, coordslinear, index):
        # index = 1 or 2
        from OpenGL import GL,GLUT
        coords = coordslinear.reshape(coordslinear.size/3, 3)
        com=np.mean(coords, axis=0)                  
        size = 0.5
        if index == 1:
            color = [0.65, 0.0, 0.0, 1.]
        else:
            color = [0.00, 0.65, 0., 1.]
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
        for i,xx in enumerate(coords):
            if i == self.ntypeA: 
                size *= 0.88 #this should be dependent on lj parameters
                if index == 1:
                    color = [0.25, 0.00, 0., 1.]
                else:
                    color = [0.00, 0.25, 0., 1.]
                GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(size,30,30)
            GL.glPopMatrix()



#
#only for testing below here
#


def run():
    #create the system object
    sys = BLJCluster(15)
    
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