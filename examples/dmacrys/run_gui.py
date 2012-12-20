import numpy as np
from pygmin.storage import savenlowest
import time
from pygmin.NEB import NEB,dimer,tstools
from pygmin.utils.rbtools import *
import dmagmin_ as GMIN
import pygmin.potentials.gminpotential as gminpot
from pygmin.takestep import generic
from pygmin.potentials.coldfusioncheck import addColdFusionCheck
import pygmin.basinhopping as bh
from pygmin.utils import dmagmin
from pygmin import defaults

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
        defaults.quenchRoutine = dmagmin.quenchCrystal
        defaults.quenchParams["tol"]=1e-4
                
    def createBasinHopping(self):
        GMIN.initialize()   
        pot = gminpot.GMINPotential(GMIN)
        coords = pot.getCoords()

        step = TestStep()
        opt = bh.BasinHopping(coords, pot, takeStep=step, quenchRoutine=dmagmin.quenchCrystal, temperature=0.4, storage=self.storage)
        addColdFusionCheck(opt)
        return opt
    
    def drawCylinder(self, X1, X2):
        from OpenGL import GL, GLU
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
        if index == 1:
            color = [[0.65, 0.0, 0.0, 1.], [0.35, 0.0, 0.0, 1.]]
        else:
            color = [[0.00, 0.65, 0., 1.], [0.00, 0.35, 0., 1.]]
            
        i=0
        for xx in coords[:-2]:
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color[i%2])
            i+=1
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
        import pygmin.utils.crystals as cr
        c1 = CoordsAdapter(nrigid=2, nlattice=6, coords=coords1)
        c2 = CoordsAdapter(nrigid=2, nlattice=6, coords=coords2)
        cr.compareStructures(c1, c2)
        
        dcom=np.mean(c2.posRigid[0], axis=0) - np.mean(c1.posRigid, axis=0)                  
        d1 = np.linalg.norm(c1.posRigid[0] - c2.posRigid[0] + dcom) + np.linalg.norm(c1.posRigid[1] - c2.posRigid[1] + dcom)
        d2 = np.linalg.norm(c1.posRigid[0] - c2.posRigid[1] + dcom) + np.linalg.norm(c1.posRigid[1] - c2.posRigid[0] + dcom)
        if(d1 > d2+1000000.):
            print "permute"
            tmp = c2.posRigid[0].copy()
            c2.posRigid[0]=c2.posRigid[1]
            c2.posRigid[1]=tmp
            tmp = c2.rotRigid[0].copy()
            c2.rotRigid[0]=c2.rotRigid[1]
            c2.rotRigid[1]=tmp
            
        for i in xrange(2):
            q1 = rot.aa2q(c1.rotRigid[i])
            q2 = rot.aa2q(c2.rotRigid[i])
            if(np.dot(q1,q2)<0):
                q2=-q2
            c2.rotRigid[i] = rot.q2aa(q2) 
            
        return coords1, coords2
    
    def createNEB(self, coords1, coords2):
        from pygmin.utils import rotations as rot
        pot = gminpot.GMINPotential(GMIN)
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
    
    def zeroEigenVecs(self, coords):
        # translational eigenvectors
        x1 = np.zeros(coords.shape)
        x2 = x1.copy()
        x3 = x1.copy()
        x1.reshape(coords.size/3,3)[0:2,0] = 1.
        x2.reshape(coords.size/3,3)[0:2,1] = 1.
        x3.reshape(coords.size/3,3)[0:2,2] = 1.

        return [x1/np.linalg.norm(x1), x2/np.linalg.norm(x2), x3/np.linalg.norm(x3)]
        
    def findTS(self, coords):
        from pygmin.optimize import fire
        pot = gminpot.GMINPotential(GMIN)
        tau = np.random.random(coords.shape) - 0.5
        tau[-6:]=0.
        from pygmin import defaults
        import pygmin.optimize.transition_state.transition_state_refinement as tsr
        defaults.quenchParams["maxstep"]=0.01
        defaults.quenchParams["tol"]=1.e-4
        defaults.quenchRoutine = fire
        ret = dimer.findTransitionState(coords+1e-2*tau, pot, direction=tau, zeroEigenVecs=self.zeroEigenVecs, tol=1e-4, maxstep=0.01)
        #ret = tsr.findTransitionState(coords+1e-2*tau, pot, tol=1e-4)
        print "TS:",ret.energy
        m1,m2 = tstools.minima_from_ts(pot.getEnergyGradient, ret.coords, ret.eigenvec, displace=1e-1)
        print "Energies: ", m1[1],ret.energy,m2[1]
        return [ret.coords,ret.energy],m1,m2
               
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(CrystalSystem)
