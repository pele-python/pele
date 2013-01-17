# utils 
import numpy as np
import tempfile
import os 

# pygmin  
from pygmin.systems import BaseSystem
from pygmin.mindist import ExactMatchAtomicCluster, MinPermDistAtomicCluster
from pygmin.transition_states import orthogopt
from pygmin.transition_states import InterpolatedPathDensity, NEB, create_NEB
from pygmin.landscape import smoothPath
from pygmin.systems import BaseParameters
from pygmin.utils.elements import elements

# OpenMM related 
from pygmin.potentials import OpenMMAmberPotential 
from simtk.unit import angstrom as openmm_angstrom

# GMIN potential  
from pygmin.potentials import GMINAmberPotential

__all__ = ["AMBERSystem_GMIN", "AMBERSystem_OpenMM"]

class AMBERBaseSystem(BaseSystem):
    """
    System class for biomolecules using AMBER ff. 
    
    Sets up using prmtop and inpcrd files used in Amber GMIN and Optim. 
    
    Potential parameters (e.g. non-bonded cut-offs are set in    
    
    TODO:   
    
    Parameters
    ----------
    prmtopFname   : str 
        prmtop file name 
    
    inpcrdFname   : str
        inpcrd file name 
        
    See Also
    --------
    BaseSystem
    """
    
    def __init__(self, prmtopFname, inpcrdFname ):
        
        super(AMBERBaseSystem, self).__init__()
        
        self.natoms = self.potential.prmtop.topology._numAtoms  
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.
        
        self.params.takestep_random_displacement = BaseParameters()
        self.params.takestep_random_displacement.stepsize = 2.
        
        self.set_params(self.params)
    
    def set_params(self, params):
        """set default parameters for the system"""
        
        #set NEBparams
        NEBparams = params.double_ended_connect.local_connect_params.NEBparams
        NEBparams.iter_density = 15.
        NEBparams.image_density = 5.
        NEBparams.max_images = 100.
        NEBparams.k = 10.
        if False: #use fire
            from pygmin.optimize import fire
            NEBparams.quenchRoutine = fire
        else: #use lbfgs
            NEBparams.NEBquenchParams.maxErise = 0.5
            NEBparams.NEBquenchParams.maxstep = .1
        NEBparams.NEBquenchParams.tol = 1e-2
            
        
        
        #set transition state search params
        tsSearchParams = params.double_ended_connect.local_connect_params.tsSearchParams
        tsSearchParams.nsteps = 200
        tsSearchParams.lowestEigenvectorQuenchParams.nsteps = 100
        tsSearchParams.lowestEigenvectorQuenchParams.tol = 0.001
        tsSearchParams.tangentSpaceQuenchParams.maxstep = .1
        tsSearchParams.nfail_max = 1000
        
        
        tsSearchParams.nsteps_tangent1 = 5
        tsSearchParams.nsteps_tangent2 = 100
        tsSearchParams.max_uphill_step = .3
        
        #control the output
        tsSearchParams.verbosity = 0
        NEBparams.NEBquenchParams.iprint = 50
        tsSearchParams.lowestEigenvectorQuenchParams.iprint = -50
        tsSearchParams.tangentSpaceQuenchParams.iprint = -5
        tsSearchParams.iprint = 10
        
#        self.params.double_ended_connect.local_connect_params.pushoff_params.verbose = True
#        self.params.double_ended_connect.local_connect_params.pushoff_params.stepmin = 1e-3
#        self.params.double_ended_connect.local_connect_params.pushoff_params.gdiff = 100.
#        #self.params.double_ended_connect.local_connect_params.pushoff_params.quenchRoutine = fire
            
    def __call__(self):
        return self 
    
    def get_potential(self):
        return self.potential 

    def get_basinhopping(self, *args, **kwargs):
        return BaseSystem.get_basinhopping(self, *args, insert_rejected=True, **kwargs)
    
    def get_takestep(self):
        from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
        
        # todo: hardcoded stepsize etc 
        takeStepRnd   = RandomDisplacement( **self.params.takestep_random_displacement )
        tsAdaptive = AdaptiveStepsizeTemperature(takeStepRnd, interval=50, verbose=False)
        return tsAdaptive     
    
    def get_random_configuration(self):
        """a starting point for basinhopping, etc."""
        from simtk.openmm.app import pdbfile as openmmpdbReader
        pdb = openmmpdbReader.PDBFile('coords.pdb')  # todo: coords.pdb is hardcoded 
        
        coords = pdb.getPositions() / openmm_angstrom
        coords = np.reshape(np.transpose(coords), 3*len(coords), 1)
        return coords 

    def get_permlist(self):
        from pygmin.utils import amberPDB_to_permList
        
        # todo: - file name coordsModTerm.pdb is hardcoded, derive from coords.pdb  
        #       - coordsModTerm.pdb should have prefix N for N terminal residue and prefix C for C terminal      
        
#        return [[0, 2, 3],    [11, 12, 13],     [19, 20, 21] ]
        if os.path.exists('coordsModTerm.pdb'):
            plist = amberPDB_to_permList.amberPDB_to_permList('coordsModTerm.pdb')
            return plist     
        else:
            print 'amberSystem: coordsModTerm.pdb not found.'    
            return []                     
        

    def get_mindist(self):
        permlist = self.get_permlist()
        
        return MinPermDistAtomicCluster(permlist=permlist, niter=10, can_invert=False)


    def createNEB(self, coords1, coords2):
        pot = self.get_potential()
        NEBparams = self.params.double_ended_connect.local_connect_params.NEBparams
        return create_NEB(pot, coords1, coords2, verbose=True, **NEBparams)

    def get_orthogonalize_to_zero_eigenvectors(self):
        return orthogopt
    
    def get_compare_exact(self, **kwargs):
        permlist = self.get_permlist()
        return ExactMatchAtomicCluster(permlist=permlist, **kwargs)

    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smoothPath(path, mindist, **kwargs)

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
        com=np.mean(coords, axis=0) 
                    
        # draw atoms as spheres      
        for i in self.potential.prmtop.topology.atoms():             
            atomElem = i.name[0]         
            atomNum  = i.index         
            x = coords[atomNum] - com 
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])            
            col = elements[atomElem]['color']
            # scaling down the radius by factor of 5, else the spheres fuse into one another 
            rad = elements[atomElem]['radius']/5  
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, col)            
            GLUT.glutSolidSphere(rad,30,30)
            GL.glPopMatrix()            
                    
        # draw bonds  
        for atomPairs in self.potential.prmtop.topology.bonds():
            # note that atom numbers in topology start at 0 
            xyz1 = coords[atomPairs[0]] - com  
            xyz2 = coords[atomPairs[1]] - com 
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
                
        #create the temporary file
        suffix = ".pdb"
        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        fname = f.name
        
        from simtk.openmm.app import pdbfile as openmmpdb

        #write the coords into pdb file
        from pygmin.mindist import CoMToOrigin
        ct = 0 
        for coords in coordslist:
            ct = ct + 1 
            coords = CoMToOrigin(coords.copy())
            self.potential.copyToLocalCoords(coords) 
#            openmmpdb.PDBFile.writeFile(self.potential.prmtop.topology , self.potential.localCoords * openmm_angstrom , file=sys.stdout, modelIndex=1)
            openmmpdb.PDBFile.writeModel(self.potential.prmtop.topology , self.potential.localCoords * openmm_angstrom , file=f, modelIndex=ct)
                        
        print "closing file"
        f.flush()
                
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

    def test_potential(self, pdbfname ):
        """ tests amber potential for pdbfname 
        
        Input
        -----
            pdbfname = full path to pdb file
             
        """                       
        # read a conformation from pdb file
        print 'reading conformation from coords.pdb' 
        from simtk.openmm.app import pdbfile as openmmpdb
        pdb = openmmpdb.PDBFile(pdbfname)        
        coords = pdb.getPositions() / openmm_angstrom   
        coords = np.reshape(np.transpose(coords), 3*len(coords), 1)
                       
        e = self.potential.getEnergy(coords)
        print 'Energy (kJ/mol) = '
        print e
        
        e, g = self.potential.getEnergyGradient(coords)
        gnum = self.potential.NumericalDerivative(coords, eps=1e-6)
    
        print 'Energy (kJ/mol) = '
        print e 
        print 'Analytic Gradient = '
        print g[1:3] 
        print 'Numerical Gradient = '
        print gnum[1:3] 
    
        print 'Num vs Analytic Gradient =' 
        print np.max(np.abs(gnum-g)), np.max(np.abs(gnum))
        print np.max(np.abs(gnum-g)) / np.max(np.abs(gnum))
        
    def test_connect(self, database):
        #connect the all minima to the lowest minimum
        minima = database.minima()
        min1 = minima[0]
        
        for min2 in minima[1:]:
            connect = self.get_double_ended_connect(min1, min2, database)
            connect.connect        
            
    def test_disconn_graph(self,database):
        from pygmin.utils.disconnectivity_graph import DisconnectivityGraph
        from pygmin.landscape import Graph
        import matplotlib.pyplot as plt
        graph = Graph(database).graph
        dg = DisconnectivityGraph(graph, nlevels=3, center_gmin=True)
        dg.calculate()
        dg.plot()
        plt.show()
            
    def test_BH(self,db):
                        
        from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
        takeStepRnd   = RandomDisplacement( stepsize=2 )
        tsAdaptive = AdaptiveStepsizeTemperature(takeStepRnd, interval=10, verbose=True)
    
        self.params.basinhopping["temperature"] = 10.0 
         
    #    sys.params.basinhopping["insert_rejected"] = True      
        # todo - how do you save N lowest?    
    
    #    bh = sys.get_basinhopping(database=database, takestep = tsAdaptive, coords=x0)    
        bh = self.get_basinhopping(database=db, takestep = takeStepRnd)     
                  
        print 'Running BH .. '
        bh.run(20)
            
        print "Number of minima found = ", len(db.minima())
        min0 = db.minima()[0]
        print "lowest minimum found has energy = ", min0.energy
    
                
    def test_mindist(self, db):
        m1, m2 = db.minima()[:2]
        mindist = sys.get_mindist()
        dist, c1, c2 = mindist(m1.coords, m2.coords)
        print "distance", dist

#  ============= Define GMIN and OpenMM specific classes as parents of AmberBaseClass 

class AMBERSystem_GMIN(AMBERBaseSystem):
    def __init__(self, prmtopFname, inpcrdFname ):
        self.potential        = GMINAmberPotential.GMINAmberPotential(prmtopFname, inpcrdFname)
        super(AMBERSystem_GMIN, self).__init__(prmtopFname, inpcrdFname)

class AMBERSystem_OpenMM(AMBERBaseSystem):
    def __init__(self, prmtopFname, inpcrdFname ):
        self.potential    = OpenMMAmberPotential.OpenMMAmberPotential(prmtopFname, inpcrdFname)
        super(AMBERSystem_OpenMM, self).__init__(prmtopFname, inpcrdFname)

#===================================================================================


if __name__ == "__main__":
    
    # create new amber system
    print '----------------------------------'
    print 'GMIN POTENTIAL' 
    sysGMIN   = AMBERSystem_GMIN('coords.prmtop', 'coords.inpcrd')        
    sysGMIN.test_potential('coords.pdb')

    print 'OPENmm POTENTIAL' 
    sysOpenMM  = AMBERSystem_OpenMM('coords.prmtop', 'coords.inpcrd')
    sysOpenMM.test_potential('coords.pdb')
    
    # load existing database 
    from pygmin.storage import Database
    dbcurr = Database(db="aladipep.db")
    
#    dbcurr.removeMinimum(1)
    
#    for i in range(1):
#        dbcurr.removeMinimum( dbcurr.getMinimum(i))
    
    print "---------id, minener"
    
    for minimum in dbcurr.minima():
        print minimum._id, minimum.energy    

    print "---------id, m1_id, m2_id, tsener"
    for ts in dbcurr.transition_states() :
        print ts._id, ts._minimum1_id, ts._minimum2_id,  ts.energy      
                
    # create new database  
    # dbcurr = sysOpenMM.create_database(db=dbcurr)    

    # connect to existing db 
#    sysOpenMM.create_database(db=dbcurr)    
    
#    for i in db.minima:
#        print i         
    
    # ------- TEST gui 
    from pygmin.gui import run as gr    
    gr.run_gui(sysOpenMM, db="aladipep.db")
    
#    # ------ Test potential 
#    sys.test_potential('coords.pdb')
#    
    # ------ BH 
#    sysOpenMM.test_BH(dbcurr)
#    # ------- Connect runs 
#    sys.test_connect(db)  
#    
#    # ------- Disconn graph  
#    sys.test_disconn_graph(db)  
#    
#    # ------- Test mindist  
#    sys.test_mindist( db)
#    





