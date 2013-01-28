import numpy as np

from pygmin.potentials import GMINPotential
import gmin_ as GMIN
from pygmin.angleaxis import RBTopology
from copy import deepcopy
from pygmin.utils import rotations
from pygmin.takestep import RotationalDisplacement
from pygmin.systems import BaseSystem, dict_copy_update, BaseParameters
from pygmin.transition_states import NEB, InterpolatedPathDensity

from pygmin.optimize import fire, mylbfgs
from pygmin import defaults

from pygmin.angleaxis.aamindist import *
from pygmin.angleaxis import MinPermDistAACluster, ExactMatchAACluster
from pygmin.angleaxis.aautils import TakestepAA
from pygmin.landscape import smoothPath
from pygmin.utils.elements import elements

class AASystem(BaseSystem):
    def __init__(self):
        BaseSystem.__init__(self)
                
        self.aasystem = self.setup_aatopology()
                
        self.params.basinhopping["temperature"]=8.
        self.params.takestep["translate"]=0.0
        self.params.takestep["rotate"]=1.6
        
        self.params.double_ended_connect.local_connect_params.nrefine_max = 5
        
        NEBparams = self.params.double_ended_connect.local_connect_params.NEBparams
        self.params.double_ended_connect.local_connect_params.NEBparams.distance=self.aasystem.neb_distance

        NEBparams.max_images=200
        NEBparams.image_density=3.0
        NEBparams.iter_density=10
        NEBparams.k = 400.
        NEBparams.adjustk_freq = 5
        NEBparams.reinterpolate = 50
        NEBparams.adaptive_nimages = True
        NEBparams.adaptive_niter = True
        NEBparams.interpolator=self.aasystem.interpolate
        NEBparams.verbose = -1
        quenchParams = NEBparams.NEBquenchParams
        #quenchParams["nsteps"] = 1000
#        quenchParams["iprint"] = -1
#        quenchParams["maxstep"] = 0.1
#        quenchParams["maxErise"] = 1000
#        quenchParams["tol"] = 1e-6
#        
        
        tsSearchParams = self.params.double_ended_connect.local_connect_params.tsSearchParams

        tsSearchParams["nfail_max"]=20
    
    def setup_aatopology(self):
        raise NotImplementedError

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
        return TakestepAA(self.aasystem, **kwargs)
    
    def get_compare_exact(self, **kwargs):
        return ExactMatchAACluster(self.aasystem, accuracy=0.1, **kwargs)
    
    def get_mindist(self, **kwargs):
        return MinPermDistAACluster(self.aasystem,accuracy=0.1, **kwargs)
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return self.aasystem.orthogopt
                
    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smoothPath(path, mindist, interpolator=self.aasystem.interpolate, **kwargs)
    
class RBSystem(AASystem):
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
        from OpenGL import GL, GLUT    
        coords = self.aasystem.to_atomistic(rbcoords)
        com=np.mean(coords, axis=0)
        self.aasystem.sites
        i=0                  
        for atom_type, xx in zip(self.atom_types, coords):
            color = [1.0, 0.0, 0.0]
            radius = 0.3
            if elements.has_key(atom_type):
                color = elements[atom_type]["color"]
                radius = elements[atom_type]["radius"]*self.render_scale
            if index == 2:
                color = [0.5, 1.0, .5]                
            
            i+=1
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(radius,30,30)
            GL.glPopMatrix()
       
        color = [1.0, 1.0, 1.0]
        if index == 2:
            color = [0.5, 1.0, .5]                
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
                
        if hasattr(self, "draw_bonds"):
            for i1, i2 in self.draw_bonds:
                self.drawCylinder(coords[i1]-com, coords[i2]-com)
