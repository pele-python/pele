"""
System class for biomolecules using AMBER ff. 
    
Set up using prmtop and inpcrd files used in Amber GMIN and Optim. 
    
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

# utils 
import numpy as np
import tempfile
import os
import shutil

# pygmin  
from pygmin.systems import BaseSystem
from pygmin.mindist import ExactMatchAtomicCluster, MinPermDistAtomicCluster
from pygmin.transition_states import orthogopt
from pygmin.transition_states import InterpolatedPathDensity, NEB, create_NEB
from pygmin.landscape import smoothPath
from pygmin.systems import BaseParameters
from pygmin.utils.elements import elements
from pygmin.systems.spawn_OPTIM import SpawnOPTIM

# OpenMM related 
import openmm_potential  
from simtk.unit import angstrom as openmm_angstrom

__all__ = ["AMBERSystem"]

class AMBERSystem(BaseSystem):

    
    def __init__(self, prmtopFname, inpcrdFname):
        
        super(AMBERSystem, self).__init__()
        
        if os.path.exists('min.in') and os.path.exists('data') :
#            print '\nFiles min.in and data found. trying to import ambgmin_ now ..' 
            try:
                import ambgmin_
                import gmin_potential
                self.potential        = gmin_potential.GMINAmberPotential(prmtopFname, inpcrdFname)
                print '\namberSystem> Using GMIN Amber potential ..'
            except ImportError:
                # using OpenMM because ambgmin_ could not be imported 
                print '\namberSystem> Using OpenMM amber potential ..'
                self.potential    = openmm_potential.OpenMMAmberPotential(prmtopFname, inpcrdFname)
        else:
            # using OpenMM because min.in and data files not found 
            print '\namberSystem> Using OpenMM amber potential ..'
            self.potential    = openmm_potential.OpenMMAmberPotential(prmtopFname, inpcrdFname)
        
        # check for openmm version
        # data structures changed between openmm4 and 5
        # crude check - todo  
        if hasattr(self.potential.prmtop.topology._bonds,'index'):
            self.OpenMMVer = 5
        else:
            self.OpenMMVer = 4                
        
        self.set_params(self.params)
        self.natoms = self.potential.prmtop.topology._numAtoms              
        
        self.params.database.accuracy = 1e-3
        self.params.basinhopping["temperature"] = 1.
        
        self.params.takestep_random_displacement = BaseParameters()
        self.params.takestep_random_displacement.stepsize = 2.
                
        self.prmtopFname = prmtopFname
        self.inpcrdFname = inpcrdFname
 
        # atom numbers of peptide bonds       
        self.populate_peptideAtomList()
        # atom numbers of CA neighbors                
        self.populate_CAneighborList() 

        self.params.basinhopping.insert_rejected = False 
        
        self.sanitycheck = True  # todo: this should be part of params and show up in GUI 
        
        if self.sanitycheck:
            self.params.basinhopping.confCheck = [self.check_cistrans_wrapper, self.check_CAchirality_wrapper]
            self.params.double_ended_connect.conf_checks = [self.check_cistrans_wrapper_kwargs, self.check_CAchirality_wrapper_kwargs]


    def set_params(self, params):
        """set default parameters for the system"""
        
        #set NEBparams
        NEBparams = params.double_ended_connect.local_connect_params.NEBparams
        NEBparams.NEBquenchParams = BaseParameters()
#        NEBquenchParams = NEBparams.NEBquenchParams

        NEBparams.iter_density = 15.
        NEBparams.image_density = 3.5
        NEBparams.max_images = 50
        NEBparams.k = 100.
        NEBparams.adjustk_freq = 5
        if False: #use fire
            from pygmin.optimize import fire
            NEBparams.quenchRoutine = fire
        else: #use lbfgs
            NEBparams.NEBquenchParams.maxErise = 100.5
            NEBparams.NEBquenchParams.maxstep = .1
        NEBparams.NEBquenchParams.tol = 1e-2
        
        NEBparams.reinterpolate = 50
        NEBparams.adaptive_niter = True
        NEBparams.adaptive_nimages = True
        NEBparams.adjustk_freq = 50                    
        
        #set transition state search params
        tsSearchParams = params.double_ended_connect.local_connect_params.tsSearchParams
        tsSearchParams.nsteps = 200
        tsSearchParams.lowestEigenvectorQuenchParams["nsteps"] = 100
        tsSearchParams.lowestEigenvectorQuenchParams["tol"] = 0.001
        tsSearchParams.tangentSpaceQuenchParams["maxstep"] = .1
        tsSearchParams.nfail_max = 1000        
        
        tsSearchParams.nsteps_tangent1 = 5
        tsSearchParams.nsteps_tangent2 = 100
        tsSearchParams.max_uphill_step = .3
        
        #control the output
        tsSearchParams.verbosity = 0
        NEBparams.NEBquenchParams.iprint = 50
        tsSearchParams.lowestEigenvectorQuenchParams["iprint"] = -50
        tsSearchParams.tangentSpaceQuenchParams["iprint"] = -5
        tsSearchParams["iprint"] = 10
        
#        self.params.double_ended_connect.local_connect_params.pushoff_params.verbose = True
#        self.params.double_ended_connect.local_connect_params.pushoff_params.stepmin = 1e-3
#        self.params.double_ended_connect.local_connect_params.pushoff_params.gdiff = 100.
#        #self.params.double_ended_connect.local_connect_params.pushoff_params.quenchRoutine = fire
            
    def __call__(self):
        return self 
    
    def get_potential(self):
        return self.potential 
        
    def get_random_configuration(self):
        """a starting point for basinhopping, etc."""
        from simtk.openmm.app import pdbfile as openmmpdbReader
        pdb = openmmpdbReader.PDBFile('coords.pdb')  # todo: coords.pdb is hardcoded 
        
        coords = pdb.getPositions() / openmm_angstrom
        coords = np.reshape(np.transpose(coords), 3*len(coords), 1)
        return coords 

    def get_metric_tensor(self, coords):
        """metric tensor for all masses m_i=1.0 """
        
        print 'amberSystem> setting up mass matrix for normal modes'
        
        # return np.identity(coords.size)
        massMatrix_tmp = np.identity(coords.size) 
        
        # get masses from 'elements' file   
        for i in self.potential.prmtop.topology.atoms():             
            atomNum  = i.index
            atomElem = i.name[0]   # assuming elements corresponding to first character of atom name       
            m        = elements[atomElem]['mass']
            massMatrix_tmp[atomNum][atomNum] = 1/m        
                    
        return massMatrix_tmp 
   
    def get_permlist(self):
        import pdb2permlist
                
        #return [[0, 2, 3], [11, 12, 13], [19, 20, 21]  ] # aladipep 
        #return [[0, 2, 3], [11, 12, 13], [21, 22, 23], [31, 32, 33], [41, 42, 43], [49,50,51]] # tetraala 
            
        if os.path.exists('coordsModTerm.pdb'):
            print '\namberSystem> constructing perm list from coordsModTerm.pdb'
            print '   (see comments in amberPDB_to_permList.py)'

            plist = pdb2permlist.pdb2permList('coordsModTerm.pdb')
            
            print '\namberSystem> Groups of permutable atoms (atom numbers start at 0) = '
            for i in plist: 
                print i
        
            return plist
        else:
            print 'amberSystem> coordsModTerm.pdb not found. permlist could not be created.'    
            return []                     
        

    def get_mindist(self):
        permlist = self.get_permlist()
        
        return MinPermDistAtomicCluster(permlist=permlist, niter=10, can_invert=False)


#    def createNEB(self, coords1, coords2):
#        pot = self.get_potential()
#        NEBparams = self.params.double_ended_connect.local_connect_params.NEBparams
#        return create_NEB(pot, coords1, coords2, verbose=True, **NEBparams)

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
            if index == 2:
                col = [0.5, 1.0, .5]                
            # scaling down the radius by factor of 5, else the spheres fuse into one another 
            rad = elements[atomElem]['radius']/5  
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, col)            
            GLUT.glutSolidSphere(rad,30,30)
            GL.glPopMatrix()            
                    
        # draw bonds  
        for atomPairs in self.potential.prmtop.topology.bonds():
            # note that atom numbers in topology start at 0
            if self.OpenMMVer == 5: 
                xyz1 = coords[atomPairs[0].index] - com  
                xyz2 = coords[atomPairs[1].index] - com
            else: 
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
        
#        #set the color according to index
#        if index == 1:
#            pymol.cmd.color("red", oname)
#        else:
#            pymol.cmd.color("blue", oname)

    def get_optim_spawner(self, coords1, coords2):
        import os
        from pygmin.config import config
        optim = config.get("exec", "AMBOPTIM")
        optim = os.path.expandvars(os.path.expanduser(optim))
        print "optim executable", optim
        return AmberSpawnOPTIM(coords1, coords2, self, OPTIM=optim, tempdir=False)


    def populate_peptideAtomList(self):
        listofC = [] 
        listofO = [] 
        listofN = [] 
        listofH = [] 
            
        for i in self.potential.prmtop.topology.atoms():
            if i.name == 'C':
                listofC.append(i.index)  
            
            if i.name == 'O':
                listofO.append(i.index)
                  
            if i.name == 'N':
                listofN.append(i.index)
                  
            if i.name == 'H':
                listofH.append(i.index)          
        
        #print listofC     
        #print listofO     
        #print listofN     
        #print listofH     
        
        # atom numbers of peptide bond 
        self.peptideBondAtoms = [] 
        
        for i in listofC: 
            if listofO.__contains__(i+1) and listofN.__contains__(i+2) and listofH.__contains__(i+3): 
                self.peptideBondAtoms.append([i,i+1,i+2,i+3]) 
        
        print '\namberSystem> Peptide bond atom numbers (C,O,N,H, in order):  '
        for i in self.peptideBondAtoms:
            print i               

    def populate_CAneighborList(self):
        listofCA = [] 
        listofC = [] 
        listofN = [] 
        listofCB = [] 
            
        for i in self.potential.prmtop.topology.atoms():
            if i.name == 'CA':
                listofCA.append(i.index)  
            
            if i.name == 'C':
                listofC.append(i.index)
                  
            if i.name == 'N':
                listofN.append(i.index)
                  
            if i.name == 'CB':
                listofCB.append(i.index)  
                
        #print '---amberSystem> list of CA = ' , listofCA              
        #print '---amberSystem> list of  C = ' , listofC              
        #print '---amberSystem> list of  N = ' , listofN              
        #print '---amberSystem> list of CB = ' , listofCB      
        
        # atom numbers of peptide bond 
        self.CAneighborList = [] 
        
        for i in listofCA:
            # find atoms bonded to CA 
            neighborlist = []                        
                 
            for b in self.potential.prmtop.topology.bonds(): 
                print b                                
                if b[0] == i:
                    neighborlist.append(b[1])
                if b[1] == i:
                    neighborlist.append(b[0])
# Commented, since this stuff doesn't seem to work at the moment...
#                if self.OpenMMVer == 5 : 
#                    # openmm5  
#                    if b[0].index == i:
#                        neighborlist.append(b[1].index)
#                    if b[1].index == i:
#                        neighborlist.append(b[0].index) 
#                else:   # openmm4 
#                    if b[0].index == i:
#                        neighborlist.append(b[1].index)
#                    if b[1].index == i:
#                        neighborlist.append(b[0].index) 
                    
            
            #print '---bonds = ', b[0].index , b[1].index 
            #print '---amberSystem> atoms bonded to CA ',i, ' = ', neighborlist    
            nn = [i] 
            # append C (=O) 
            for n in neighborlist: 
                if listofC.__contains__(n):
                    nn.append(n) 
        
            # append CB  
            for n in neighborlist: 
                if listofCB.__contains__(n):
                    nn.append(n) 
        
            # append N  
            for n in neighborlist: 
                if listofN.__contains__(n):
                    nn.append(n) 
        
            self.CAneighborList.append(nn) 

        # atoms numbers start at 0             
        print '\namberSystem> CA neighbors atom numbers (CA,C(=O),CB, N, in order):  '
        for i in self.CAneighborList:
            print i               

    def check_cistrans_wrapper_kwargs(self, coords=None, **kwargs):
        print 'in check_cistrans_wrapper_kwargs'
        return self.check_cistrans(coords)

    def check_cistrans_wrapper(self, energy, coords, **kwargs):
        return self.check_cistrans(coords)

    def check_cistrans(self, coords):
        """ 
        Sanity check on the isomer state of peptide bonds   
        
        Returns False if the check fails i.e. if any of the peptide bond is CIS         
        
        """
        
        import measure  
        m = measure.Measure() 
        
        isTrans = True 
        
        for i in self.peptideBondAtoms:                            
            atNum = i[0] 
            rC    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[1] 
            rO    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[2] 
            rN    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[3] 
            rH    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])
            
            # compute O-C-N-H torsion angle 
            rad, deg = m.torsion(rO,rC,rN,rH)
                        
            # print 'peptide torsion (deg) ', i, ' = ', deg 
            # check cis 
            if deg < 90 or deg > 270: 
                isTrans = False  
#                print 'CIS peptide bond between atoms ', i, ' torsion (deg) = ', deg 
                
        return isTrans  
            

    def check_CAchirality_wrapper_kwargs(self, coords=None, **kwargs):
        return self.check_CAchirality(coords)
    
    def check_CAchirality_wrapper(self, energy, coords, **kwargs):
        return self.check_CAchirality(coords)

    def check_CAchirality(self, coords):
        """ 
        Sanity check on the CA to check if it is L of D    
        
        Returns False if the check fails i.e. if any D-amino acid is present          
        
        """
        
        # print 'in check CA chirality'
        import measure  
        m = measure.Measure() 
        
        isL = True 
        
        for i in self.CAneighborList:                            
            atNum = i[0] 
            rCA   = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[1] 
            rC    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[2] 
            rCB   = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])                       
            atNum = i[3] 
            rN    = np.array( [ coords[3*atNum] , coords[3*atNum+1] , coords[3*atNum+2] ])
            
            # compute improper torsion angle between C-CA-CB and CA-CB-N 
            rad, deg = m.torsion(rC,rCA,rCB,rN)                        
            
            # check cis 
            if deg < 180 :
                # this condition was found by inspection of structures todo   
                isL = False  
#                print 'chiral state of CA atom ', i[0], ' is D' 
#                print 'CA improper torsion (deg) ', i, ' = ', deg 

        return isL  


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
            
    def test_BH(self,db,nsteps):
                        
        from pygmin.takestep import RandomDisplacement, AdaptiveStepsizeTemperature
        takeStepRnd   = RandomDisplacement( stepsize=2 )
        tsAdaptive = AdaptiveStepsizeTemperature(takeStepRnd, interval=10, verbose=True)
    
        self.params.basinhopping["temperature"] = 10.0 
         
        # todo - how do you save N lowest?    
    
        bh = self.get_basinhopping(database=db, takestep = takeStepRnd)     
        bh = self.get_basinhopping(database=db, takestep = tsAdaptive)     
                  
        print 'Running BH .. '
        bh.run(nsteps)
            
        print "Number of minima found = ", len(db.minima())
        min0 = db.minima()[0]
        print "lowest minimum found has energy = ", min0.energy
    
                
    def test_mindist(self, db):
        m1, m2 = db.minima()[:2]
        mindist = sys.get_mindist()
        dist, c1, c2 = mindist(m1.coords, m2.coords)
        print "distance", dist

class AmberSpawnOPTIM(SpawnOPTIM):
    def __init__(self, coords1, coords2, sys, **kwargs):
        super(AmberSpawnOPTIM, self).__init__(coords1, coords2, **kwargs)
        self.sys = sys
    
    def write_odata_coords(self, coords, fout):
        pass

    def write_perm_allow(self, fname):
        permallow = self.make_permallow_from_permlist(self.sys.get_permlist())
        with open(fname, "w") as fout:
            fout.write(permallow)
    
    def write_additional_input_files(self, rundir, coords1, coords2):
        #write start
        with open(rundir + "/start", "w") as fout:
            for xyz in coords1.reshape(-1,3):
                fout.write( "%f %f %f\n" % tuple(xyz))
        
        #write coords.prmtop and coords.inpcrd
        shutil.copyfile(self.sys.prmtopFname, rundir + "/coords.prmtop")
        shutil.copyfile(self.sys.inpcrdFname, rundir + "/coords.inpcrd")
        min_in = """
STOP
 &cntrl
  imin   = 1,
  ncyc = 1,
  maxcyc = 1,
  igb = 0,
  ntb    = 0,
  cut    = 999.99,
  rgbmax = 25.0,
  ifswitch = 1
 /
"""
        with open(rundir + "/min.in", "w") as fout:
            fout.write(min_in)
            
    
    def write_odata(self, fout):
        odatastr = """
DUMPALLPATHS

UPDATES 6000
NEWCONNECT 15 3 2.0 20.0 30 0.5
CHECKCHIRALITY
comment PATH dumps intermediate conformations along the path
PATH 100 1.0D-2
COMMENT NEWNEB 30 500 0.01
NEBK 10.0
comment DUMPNEBXYZ
AMBERIC
comment AMBERSTEP
DIJKSTRA EXP
DUMPALLPATHS
REOPTIMISEENDPOINTS
COMMENT MAXTSENERGY -4770.0
EDIFFTOL  1.0D-4
MAXERISE 1.0D-4 1.0D0
GEOMDIFFTOL  0.05D0
BFGSTS 500 10 100 0.01 100
NOIT
BFGSMIN 1.0D-6
PERMDIST
MAXSTEP  0.1
TRAD     0.2
MAXMAX   0.3
BFGSCONV 1.0D-6
PUSHOFF 0.1
STEPS 800
BFGSSTEPS 2000
MAXBFGS 0.1
NAB start
"""
        fout.write(odatastr)
        fout.write("\n")

# ============================ MAIN ================================ 

if __name__ == "__main__":
    
    # create new amber system    
    sysAmb  = AMBERSystem('/home/ss2029/WORK/PyGMIN/examples/amber/coords.prmtop', '/home/ss2029/WORK/PyGMIN/examples/amber/coords.inpcrd')
    
    # load existing database 
    from pygmin.storage import Database
    dbcurr = Database(db="/home/ss2029/WORK/PyGMIN/examples/amber/aladipep.db")
                        
    
    coords = sysAmb.get_random_configuration()
    aa = sysAmb.get_metric_tensor(coords)
    
    exit()  
    # ------- TEST gui 
    from pygmin.gui import run as gr    
    gr.run_gui(sysAmb)
    
    # ------ Test potential 
    sysAmb.test_potential('coords.pdb')
    
    # ------ BH 
    sysAmb.test_BH(dbcurr)
    
    # ------- Connect runs 
    sysAmb.test_connect(dbcurr)  
    
    # ------- Disconn graph  
    sysAmb.test_disconn_graph(dbcurr)  
    
    # ------- Test mindist  
    sysAmb.test_mindist( dbcurr)
    
