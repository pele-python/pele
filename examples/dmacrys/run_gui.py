import numpy as np
from pygmin.storage import savenlowest
import time
from pygmin.NEB import NEB
from pygmin.utils.rbtools import *
import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
from pygmin.optimize import quench
from pygmin.takestep import generic
from pygmin.potentials.coldfusioncheck import addColdFusionCheck
import pygmin.basinhopping as bh
from pygmin.utils import dmagmin

class TestStep(generic.TakestepInterface):
    def takeStep(self, coords, **kwargs):
        from pygmin.takestep import buildingblocks as bb
        print coords.size
        ca = CoordsAdapter(nrigid=2, nlattice=6, coords=coords)
        bb.rotate(1.6, ca.rotRigid)
        ca.lattice*=1.2
        
class CrystalSystem:
    def __init__(self):
        self.storage = savenlowest.SaveN(10)
        GMIN.initialize()
        
    def createBasinHopping(self):
        GMIN.initialize()   
        pot = gminpot.GMINPotental(GMIN)
        coords = pot.getCoords()

        step = TestStep()
        opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=dmagmin.quenchCrystal, temperature=0.4, storage=self.storage)
        addColdFusionCheck(opt)
        return opt
    
    def drawCylinder(self, X1, X2):
        from OpenGL import GL,GLUT, GLU
        z = np.array([0.,0.,1.]) #default cylinder orientation
        p = X2-X1 #desired cylinder orientation
        r = np.linalg.norm(p)
        t = np.cross(z,p)  #angle about which to rotate
        a = np.arccos( np.dot( z,p) / r ) #rotation angle
        a *= (180. / np.pi)  #change units to angles
        if(np.dot(t,t) < 1e-6):
             t = [1.,0.,0.]
        GL.glPushMatrix()
        GL.glTranslate( X1[0], X1[1], X1[2] )
        if (np.linalg.norm(t) > 1e-3):
            GL.glRotate( a, t[0], t[1], t[2] )
        g=GLU.gluNewQuadric()
        GLU.gluCylinder(g, .1,0.1,r,30,30)  #I can't seem to draw a cylinder
        GL.glPopMatrix()

    def draw(self, coords_rigid, index):
        from OpenGL import GL,GLUT
        coords=np.zeros([GMIN.getNAtoms(), 3])
        GMIN.toAtomistic(coords.reshape(coords.size), coords_rigid)
        #coords = coords.reshape(GMIN.getNAtoms, 3)
        com=np.mean(coords, axis=0)                  
        for xx in coords[:-2]:
            x = xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.5,30,30)
            GL.glPopMatrix()
        
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, [0.,0.,1.,1.])
        from pygmin.utils import lattice
        l = lattice.lowerTriangular(coords_rigid[-6:])
        mb = l[:,0] + l[:,1] + l[:,2]
        GL.glPushMatrix()            
        GL.glTranslate(-0.5*mb[0],-0.5*mb[1],-0.5*mb[2])
        
        for i in xrange(3):
            self.drawCylinder([0.,0.,0.], l[:,i])
        for i in xrange(1,3):
            self.drawCylinder(l[:,0], l[:,0] + l[:,i])
        mb = l[:,0] + l[:,1] + l[:,2]
        for i in xrange(0,3):
            self.drawCylinder(mb, mb - l[:,i])
        for i in xrange(1,3):
            self.drawCylinder(mb-l[:,0], mb - l[:,i]-l[:,0])
        self.drawCylinder(l[:,1], l[:,1] + l[:,0])            
        self.drawCylinder(l[:,2], l[:,2] + l[:,0])            
        GL.glPopMatrix()
        
    def Align(self, coords1, coords2):
        import pygmin.utils.rotations as rot
        c1 = CoordsAdapter(nrigid=2, nlattice=6, coords=coords1)
        c2 = CoordsAdapter(nrigid=2, nlattice=6, coords=coords2)
        for i in xrange(2):
            q1 = rot.aa2q(c1.rotRigid[i])
            q2 = rot.aa2q(c2.rotRigid[i])
            if(np.dot(q1,q2)<0):
                q2=-q2
            c2.rotRigid[i] = rot.q2aa(q2) 
            
        return coords1, coords2
    
    def createNEB(self, coords1, coords2):
        from pygmin.utils import rotations as rot
        pot = gminpot.GMINPotental(GMIN)
        neb = NEB.NEB(coords1, coords2, pot, k = 100. ,nimages=20)
        # replace rotational interpolation by slerp
        c1 = CoordsAdapter(nrigid=2, nlattice=6, coords=neb.coords[0, :])
        c2 = CoordsAdapter(nrigid=2, nlattice=6, coords=neb.coords[-1, :])

        for i in xrange(1, neb.nimages):
            ci = CoordsAdapter(nrigid=2, nlattice=6, coords=neb.coords[i, :])
            t = float(i) / float(neb.nimages)
            for j in xrange(2):
                ci.rotRigid[j] =  rot.q2aa(rot.q_slerp(
                                       rot.aa2q(c1.rotRigid[j]),
                                       rot.aa2q(c2.rotRigid[j]), t))
        return neb
               
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(CrystalSystem)
