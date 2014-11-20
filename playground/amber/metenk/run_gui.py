import numpy as np

import pele.utils.readAmberParam as readAmb
import ambgmin_ as GMIN
import pele.potentials.gminpotential as gminpot
import pele.basinhopping as bh
from pele.storage import savenlowest
from pele.NEB import NEB
from pele.takestep import displace


class molSystem:
    def __init__(self):
        self.storage = savenlowest.SaveN(10)
        GMIN.initialize()
#        self.bondList = bondList 
        
    def createBasinHopping(self):
        GMIN.initialize()   
        pot = gminpot.GMINPotental(GMIN)
        coords = pot.getCoords()

        step = displace.RandomDisplacement()
        opt = bh.BasinHopping(coords, pot, takeStep=step, temperature=0.4, storage=self.storage)
        return opt
    
    def drawCylinder(self, X1, X2):
        from OpenGL import GL, GLU
        z = np.array([0.,0.,1.]) #default cylinder orientation
        p = X2-X1 #desired cylinder orientation
        r = np.linalg.norm(p)
        t = np.cross(z,p)  #angle about which to rotate
        a = np.arccos( np.dot( z,p) / r ) #rotation angle
        a *= (180. / np.pi)  #change units to angles
        GL.glPushMatrix()
        GL.glTranslate( X1[0], X1[1], X1[2] )
        GL.glRotate( a, t[0], t[1], t[2] )
        g=GLU.gluNewQuadric()
        GLU.gluCylinder(g, .1,0.1,r,30,30)  #I can't seem to draw a cylinder
        GL.glPopMatrix()
        
    def draw(self, coordsl, index):
        from OpenGL import GL,GLUT
        coords=coordsl.reshape(coordsl.size/3,3)
        #coords = coords.reshape(GMIN.getNAtoms, 3)
        com=np.mean(coords, axis=0)                  
        for xx in coords:
            x = xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.3,30,30)
            GL.glPopMatrix()

        # get bond list from amber params 
        mol = readAmb.readAmberParam()
        mol.populateBondConn() 
        
        # draw bonds  
        for atomPairs in mol.bondConn:
            xyz1 = coords[atomPairs[0]-1] - com  
            xyz2 = coords[atomPairs[1]-1] - com 
            self.drawCylinder(xyz1, xyz2)            
                     
    def createNEB(self, coords1, coords2):
        pot = gminpot.GMINPotental(GMIN)
        return NEB.NEB(coords1, coords2, pot, k = 100. ,nimages=20)

               
if __name__ == "__main__":
    import pele.gui.run as gr
    gr.run_gui(molSystem)
