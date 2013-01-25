import numpy as np
from pygmin import gui
import oxdnagmin_ as GMIN
import pygmin.gui.run as gr
from pygmin.utils.rbtools import CoordsAdapter
from pygmin import takestep
from pygmin.utils import rotations
from pygmin.potentials import GMINPotential
import pygmin.basinhopping as bh
from pygmin.transition_states import NEB, InterpolatedPath
        
# This is the takestep routine for OXDNA. It is a standard rigid body takestep
# routine, but I put it here to be able to start modifying it
class OXDNATakestep(takestep.TakestepInterface):
    def __init__(self, displace=0.5, rotate=0.5):
        self.displace = displace
        self.rotate = rotate
        
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        ca.posRigid[:] += 2.*self.displace*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        takestep.rotate(self.rotate, ca.rotRigid)
        
    # this is necessary for adaptive step taking
    def scale(self, factor):
        self.rotate *= factor
        self.displace *= factor

    @property
    def stepsize(self):
        return [self.rotate, self.displace]

# this class should generate a fully random configuration
class OXDNAReseed(takestep.TakestepInterface):
    def __init__(self, radius=3.0):
        self.radius = radius
    
    def takeStep(self, coords, **kwargs):
        # easy access to coordinates
        ca = CoordsAdapter(nrigid=coords.size/6, coords = coords)
        
        # random displacement for positions
        ca.posRigid[:] = 2.*self.radius*(np.random.random(ca.posRigid.shape)-0.5)
        
        # random rotation for angle-axis vectors
        for rot in ca.rotRigid:
            rot[:] = rotations.random_aa()


class OXDNASystem(gui.GUISystem):
    def create_basinhopping(self):
        potential = GMINPotential(GMIN)
        coords = potential.getCoords()
        # genearte a random configuration
        OXDNAReseed().takeStep(coords)
        
        opt = bh.BasinHopping(coords,potential,
                          temperature=0.1, takeStep=self.create_takestep(),
                          outstream = None)
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
        GLU.gluCylinder(g, .05,0.05,r,30,30)  #I can't seem to draw a cylinder
        GL.glPopMatrix()


    def draw(self, coordslinear, index):
        ca = CoordsAdapter(nrigid=coordslinear.size/6,coords=coordslinear)
        from OpenGL import GL,GLUT
        com=np.mean(ca.posRigid, axis=0)                  
        for xx, rot in zip(ca.posRigid, ca.rotRigid):
            p = np.dot(rotations.aa2mx(rot), np.array([1., 0., 0.]))
            x=xx-com
            x = x - 0.4*p
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(0.3,30,30)
            p=0.8*p
            self.drawCylinder([0.,0.,0.], p)
            GL.glTranslate(p[0],p[1],p[2])
            GLUT.glutSolidSphere(0.1,30,30)
            GL.glPopMatrix()
        
    def create_takestep(self):
        group = takestep.BlockMoves()

        step1 = takestep.AdaptiveStepsize(OXDNATakestep(displace=0.0, rotate=0.5), frequency=50)
        #step2 = takestep.AdaptiveStepsize(OXDNATakestep(displace=0., rotate=0.5), frequency=50)
        group.addBlock(100, step1)
        #group.addBlock(100, step2)

        # with a generate random configuration
        genrandom = OXDNAReseed()
        # in a reseeding takestep procedure
        reseed = takestep.Reseeding(group, genrandom, maxnoimprove=2000)
        return reseed
    
    def Align(self, coords1, coords2):
        from pygmin.mindist import rmsfit,aamindist
        potential = GMINPotential(GMIN)
        
        ca1 = CoordsAdapter(nrigid=coords1.size/6, coords = coords1)
        ca2 = CoordsAdapter(nrigid=coords2.size/6, coords = coords2)
        rot = rmsfit.findrotation_kabsch(ca2.posRigid, ca1.posRigid)
        ca2.posRigid[:] = np.dot(rot,ca2.posRigid.transpose()).transpose()
        for p in ca2.rotRigid:
            p[:] = rotations.mx2aa((np.dot(rot, rotations.aa2mx(p))))
            
        print "before"
        print potential.getEnergy(coords1), potential.getEnergy(coords2)
        print ca2.rotRigid - ca1.rotRigid
        for p1,p2 in zip(ca1.rotRigid, ca2.rotRigid):
            p1[:],p2[:] = aamindist.aadistance(p1, p2)
        print "after"
        print potential.getEnergy(coords1), potential.getEnergy(coords2)
        print ca2.rotRigid - ca1.rotRigid
        path = InterpolatedPath(coords1, coords2, 40)
        print potential.getEnergy(path[0])
        print "Interpolated energies"
        for x in InterpolatedPath(coords1, coords2, 40):
            print potential.getEnergy(x)
        print "done"
        return coords1, coords2
    
    def createNEB(self, coords1, coords2):
        return NEB(InterpolatedPath(coords1, coords2, 40), GMINPotential(GMIN), k = 100.)
    
        
if __name__ == "__main__":
    GMIN.initialize()
    gr.run_gui(OXDNASystem)