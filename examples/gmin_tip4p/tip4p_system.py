import numpy as np

from pygmin.potentials import GMINPotential
import gmin_ as GMIN
from pygmin.angleaxis import RBSystem
from copy import deepcopy
import tip4p
from pygmin.utils import rotations
from pygmin.takestep import RotationalDisplacement
from pygmin.systems import BaseSystem, dict_copy_update, BaseParameters
from pygmin.transition_states import NEB, InterpolatedPathDensity

from pygmin.optimize import fire, mylbfgs
from pygmin import defaults

from pygmin.angleaxis.aamindist import *
from pygmin.mindist import MinPermDistCluster

class TIP4PSystem(BaseSystem):
    def __init__(self):
        BaseSystem.__init__(self)
        neb_params = BaseParameters()
        self.params["neb"]=neb_params
        
        neb_params["nimages"] = 15
        neb_params["k"] = 10.
        neb_params["aadist"] = True
        
        
        defaults.NEBquenchParams["nsteps"] = 200
        defaults.NEBquenchParams["iprint"] = -1
        defaults.NEBquenchParams["maxstep"] = 0.1
        #defaults.NEBquenchParams["maxErise"] = 0.1
        defaults.NEBquenchParams["tol"] = 1e-6
        defaults.NEBquenchRoutine = fire     
        
        GMIN.initialize()
        pot = GMINPotential(GMIN)
        coords = pot.getCoords()

        nrigid = coords.size / 6
        print "I have %d water molecules in the system"%nrigid

        water = tip4p.water()
        system = RBSystem()
        system.add_sites([deepcopy(water) for i in xrange(nrigid)])
        
        self.aasystem = system
        self.potential = pot
        self.nrigid = nrigid

    def get_potential(self):
        return self.potential
    
    def get_random_configuration(self):
        coords = 5.*np.random.random(6*self.nrigid)
        ca = self.aasystem.coords_adapter(coords)
        for p in ca.rotRigid:
            p = rotations.random_aa()
        return coords
    
    def get_takestep(self, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        
        default is random displacement with adaptive step size 
        adaptive temperature
        
        See Also
        --------
        pygmin.takestep
        """
        kwargs = dict_copy_update(self.params["takestep"], kwargs)                
        return RotationalDisplacement(**kwargs)
    
    def get_mindist(self, **kwargs):
        transform = TransformAngleAxisCluster(self.aasystem)
        measure = MeasureAngleAxisCluster(self.aasystem, transform)
        return MinPermDistCluster(transform = transform, measure=measure)
    
    def createNEB(self, coords1, coords2):
        pot = self.get_potential()
        #dist = np.linalg.norm(coords1- coords2)
        #if dist < 1.: dist = 1
        #image_density = 15.
        
        #path = InterpolatedPathDensity(coords1, coords2, 
        #                               distance=dist, density=image_density)
        path = tip4p.get_path(self.aasystem, coords1, coords2, self.params.neb.nimages)
                              
        if(self.params.neb.aadist):
            return NEB(path, pot, k = self.params.neb.k, distance=self.aasystem.neb_distance)
        else:
            return NEB(path, pot, k = self.params.neb.k)
    
    def drawCylinder(self, X1, X2):
        from OpenGL import GL,GLUT, GLU
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
        
    def draw(self, rbcoords, index):
        from OpenGL import GL,GLUT
        coords = self.aasystem.to_atomistic(rbcoords)
        com=np.mean(coords, axis=0)
        i=0                  
        for xx in coords:
            if(i%3 == 0):
                color = [1.0, 0.0, 0.0]
                radius = 0.35
            else:
                color = [1.0, 1.0, 1.0]
                radius = 0.3
            if index == 2:
                color = [0.5, 1.0, .5]                
            
            i+=1
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(radius,30,30)
            GL.glPopMatrix()
        
        for i in xrange(self.nrigid):
            color = [1.0, 0.0, 0.0]
            if index == 2:
                color = [0.5, 1.0, .5]                
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            self.drawCylinder(coords[3*i]-com, coords[3*i+1]-com)
            self.drawCylinder(coords[3*i]-com, coords[3*i+2]-com)
            
if __name__ == "__main__":
    import pygmin.gui.run as gr
    gr.run_gui(TIP4PSystem)