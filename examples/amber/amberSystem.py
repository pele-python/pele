import numpy as np

import tempfile

from pygmin.systems import BaseSystem
import openMMamberPot as amb  
import pygmin.utils.amber as readAmb

from simtk.unit import  angstrom 

from pygmin.mindist import minPermDistStochastic, MinDistWrapper, ExactMatchCluster
from pygmin.transition_states import orthogopt
from pygmin.transition_states import InterpolatedPathDensity, NEB

from simtk.openmm.app import pdbfile as openmmpdb
from simtk.unit import  * 

__all__ = ["AMBERsystem"]

class AMBERsystem(BaseSystem):
    """
    define the System class for a AMBER biosystem 
    
    Parameters
    ----------
    prmtopFname   = prmtop file name 
    
    inpcrdFname   = inpcrd file name 
        
    See Also
    --------
    BaseSystem
    """
    
    def __init__(self, prmtopFname, inpcrdFname ):
        
        super(AMBERsystem, self).__init__()
            
        # set-up potential                 
        self.pot    = amb.amberPot(prmtopFname, inpcrdFname) 
            
        self.natoms = self.pot.prmtop.topology._numAtoms  
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.0
        self.params.double_ended_connect.local_connect_params.NEBparams.image_density = 1.
            
    def __call__(self):
        return self 
    
    def get_potential(self):
        return self.pot 

    def get_basinhopping(self, *args, **kwargs):        
        return BaseSystem.get_basinhopping(self, *args, insert_rejected=True, **kwargs)
    
    def get_takestep(self):
        from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
        takeStepRnd   = RandomDisplacement( stepsize=2 )
        tsAdaptive = AdaptiveStepsizeTemperature(takeStepRnd, interval=10, verbose=True)
        return tsAdaptive     
    
    def get_random_configuration(self):
        """a starting point for basinhopping, etc."""
        pdb = openmmpdb.PDBFile('coords.pdb')
        
        coords = pdb.getPositions() / angstrom
        coords = np.reshape(np.transpose(coords), 3*len(coords), 1)
        return coords 

    def get_mindist(self):
        # todo - permutable hydrogens for ala dipep 
        permlist = [[0, 2, 3],    [11, 12, 13],     [19, 20, 21] ]
#        permlist = []
        
        return MinDistWrapper(minPermDistStochastic, permlist=permlist, niter=10)

    def createNEB(self, coords1, coords2):
        pot = self.get_potential()
        dist = np.linalg.norm(coords1- coords2)
        if dist < 1.: dist = 1
        try :
            image_density = self.params.double_ended_connect.local_connect_params.NEB_image_density
        except:
            image_density = 1.
            print "using image_density", image_density
        
        path = InterpolatedPathDensity(coords1, coords2, 
                                       distance=dist, density=image_density)
        return NEB(path, pot)

    def get_orthogonalize_to_zero_eigenvectors(self):
        return orthogopt
    
    def get_compare_exact(self, **kwargs):
        permlist = [range(self.natoms)]
        return ExactMatchCluster(permlist=permlist, **kwargs)
    
    #
    # ========= below here is stuff only for the gui =====================
    #
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

        
    def load_coords_pymol(self, coordslist, oname, index=1):
        """load the coords into pymol
        
        the new object must be named oname so we can manipulate it later
                        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object must be named oname so it can be manipulated
            later
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2
        
        Notes
        -----
        the implementation here is a bit hacky.  we create a temporary xyz file from coords
        and load the molecule in pymol from this file.  
        """
        #pymol is imported here so you can do, e.g. basinhopping without installing pymol
        import pymol 
        
        
        # todo -- work with 
#        #create the temporary file
#        suffix = ".pdb"
#        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
#        fname = f.name

        fname = 'temp.pdb'
        f = open(fname,'w')
        
        #write the coords into pdb file
        from pygmin.mindist import CoMToOrigin
        ct = 0 
        for coords in coordslist:
            ct = ct + 1 
            coords = CoMToOrigin(coords.copy())
            self.pot.copyToLocalCoords(coords) 
#            openmmpdb.PDBFile.writeFile(self.pot.prmtop.topology , self.pot.localCoords * angstrom , file=sys.stdout, modelIndex=1)
            openmmpdb.PDBFile.writeModel(self.pot.prmtop.topology , self.pot.localCoords * angstrom , file=f, modelIndex=ct)
                        
#        f.flush()
        print "closing file"
        f.close() 
                
        #load the molecule from the temporary file
        pymol.cmd.load(fname)
        
        #get name of the object just created and change it to oname
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        pymol.cmd.set_name(objectname, oname)
        
        #set the representation
        pymol.cmd.hide("everything", oname)
        pymol.cmd.show("lines", oname)
        
        #set the color according to index
        if index == 1:
            pymol.cmd.color("red", oname)
        else:
            pymol.cmd.color("blue", oname)

#
# ===================================================================== 
#    
if __name__ == "__main__":
    
    # create new amber system 
    sys   = AMBERsystem('coords.prmtop', 'coords.inpcrd')        
    # create new database  
    db = sys.create_database()
    
#    sys.draw()
    #create a database
    db = sys.create_database()
    
    #do a short basinhopping run
#    bh = sys.get_basinhopping(database=db, outstream=None)
#    while len(db.minima()) < 2:
#        bh.run(100)
    
    #try to connect the lowest two minima
#    min1, min2 = db.minima()[:2]
#    connect = sys.get_double_ended_connect(min1, min2, db)
#    connect.connect()

#    run()
