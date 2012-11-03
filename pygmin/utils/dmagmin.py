import dmagmin_ as GMIN
from pygmin.utils import rotations
from pygmin import takestep
from pygmin.optimize import quench
import crystals
import lattice
import numpy as np
from rbtools import CoordsAdapter
from pygmin.application import AppClusterBH
from pygmin.application.basinhopping_app import AppBasinHopping
from pygmin.potentials.coldfusioncheck import addColdFusionCheck
from pygmin.potentials.gminpotential import GMINPotential 

__all__ = ["compareMinima", "GenRandomCrystal", "quenchCrystal",
           "TakestepDMAGMIN", "AppDMAGMINBH", "DMACRYSPotential"]

class DMACRYSPotential(GMINPotential):
    ''' wrapper for dmacrys potential
    
        This is a wrapper for the DMACRYS potential in GMIN. It has some
        additional features like selecting energy units or freezing
        the lattice        
    '''  
    
    ENERGY_IN_EV = 1.0
    ENERGY_IN_KJMOL = 96.4853365
        
    def __init__(self):
        GMINPotential.__init__(self, GMIN)
        self.energy_unit = self.ENERGY_IN_KJMOL
        self.freeze_lattice = False
    
    def getEnergy(self, coords):
        return self.energy_unit * GMINPotential.getEnergy(self, coords)

    def getEnergyGradient(self, coords):        
        E, g = GMINPotential.getEnergyGradient(self, coords)
        if self.freeze_lattice:
            g[-6:]=0.
        g*=self.energy_unit
        return self.energy_unit * E, g
        
# compare 2 minima
def compareMinima(min1, min2):
    from pygmin.utils import rbtools
    ca1 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min1.coords)
    ca2 = rbtools.CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=min2.coords)
    match = crystals.compareStructures(ca1, ca2)
    return match

# generate a random crystal structure
class GenRandomCrystal(takestep.TakestepInterface):
    ''' takestep class to generate a random crystal
    
        GenRandomCrystal is a takestep class which generates a random crystal structure.
        It can be either used as a standard takestep routine to perform a random search
        or as a reseeding routine in combination with pygmin.takestep.group.Reseeding
    '''
    
    def __init__(self, coordsadapter, volume=None, shear=2., expand=4.0, expand_current=1.2, overlap=None):
        '''
        :param volume: volume bounds of the generated cell [min, max]. None to use 2* current volume
        :param shear: maximum off diagonal matrix element for lattice matrix
        :param expane: maxumum assymmetry of the cell, 0 means generate a cubic cell
        
        ''' 
        
        self.volume = volume
        self.shear = shear
        self.expand = expand
        self.coordsadapter = coordsadapter
        self.expand_current = expand_current
        self.overlap = overlap
        
    def takeStep(self, coords, **kwargs):
        ''' takeStep routine to generate random cell '''        
        ca = self.coordsadapter        
        ca.updateCoords(coords)
        
        atomistic = np.zeros(3*GMIN.getNAtoms())
        valid_configuration = False
        for i in xrange(50):
            volumeTarget = self.expand_current*lattice.volume(ca.lattice)
             
            # random box
            ca.lattice[[0,3,5]] = 1.0 + self.expand * np.random.random(3)  
            ca.lattice[[1,2,4]] = self.shear * np.random.random(3)
            
            if(self.volume != None):
                volumeTarget = self.volume[0] + (self.volume[1] - self.volume[0]) * np.random.random()
                        
            vol = lattice.volume(ca.lattice)
            ca.lattice[:] = ca.lattice * (volumeTarget / vol)**(1.0/3.0)
            GMIN.reduceCell(coords)
            
            for i in xrange(50):# first choose random positions and rotations
                for i in xrange(GMIN.getNRigidBody()):
                    ca.posRigid[i] = np.random.random()
                    ca.rotRigid[i] = rotations.random_aa()
        
                if self.overlap is None:
                    return
            
            
                GMIN.toAtomistic(atomistic, coords)
                if not crystals.has_overlap(atomistic, self.overlap):
                    return
                            
            print "Could generate valid configuration for current box, choose new box"
        raise Exception("GenRandomCrystal: failed to generate a non-overlapping configuration")
            
# special quencher for crystals
def quenchCrystal(coords, pot, **kwargs):
    ''' Special quench routine for crystals which makes sure that the final structure is a reduced cell '''
    coords, E, rms, calls = quench.lbfgs_py(coords, pot, **kwargs)
    #while(GMIN.reduceCell(coords)):
    if(GMIN.reduceCell(coords)):
        #print "Reduced cell, redo minimization"
        coords, E, rms, callsn = quench.lbfgs_py(coords, pot, **kwargs)
        calls+=callsn
    return coords, E, rms, calls

class TakestepDMAGMIN(takestep.TakestepInterface):
    def __init__(self, expand=1.0, rotate=1.6, translate=0., nmols=None, overlap_cutoff=None, max_volume=None):
        self.expand = expand
        self.rotate = rotate
        self.translate = translate
        self.nmols=nmols
        self.overlap_cutoff = overlap_cutoff
        self.max_volume = max_volume
        
    def takeStep(self, coords, **kwargs):
        from pygmin.takestep import buildingblocks as bb
        
        ca = CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=coords)
        
        #tmp = np.loadtxt("aligned.dat")
        #coords[-6:]=tmp[-6:]        
        #coords[:]=tmp[:]
        #return        
        conf_ok = False
        if(not self.overlap_cutoff is None):
            safe_coords = coords.copy()
        
        while conf_ok is False:        
            indices=None
            if(self.nmols):
                indices=[]
                indices = [np.random.randint(0,GMIN.getNRigidBody()) for i in xrange(self.nmols)]
            
            if(self.rotate != 0.):
                bb.rotate(self.rotate, ca.rotRigid, indices)
            if(self.translate != 0.):
                bb.uniform_displace(self.translate, ca.rotRigid, indices)
                
            #from pygmin.utils import lattice
            #bb.reduced_coordinates_displace(0.0, lattice.lowerTriangular(ca.lattice), ca.posRigid)
            ca.lattice[:]*=self.expand

            if(self.max_volume):
                vol = lattice.volume(ca.lattice)
                if(vol > self.max_volume):
                    ca.lattice[:]*=(self.max_volume / vol)**(1./3.)
            
            conf_ok = True
            if(self.overlap_cutoff):
                atomistic = np.zeros(3*GMIN.getNAtoms())
                GMIN.toAtomistic(atomistic, coords)                
                overlap =  crystals.has_overlap(atomistic, self.overlap_cutoff)
                if(overlap is True):
                    conf_ok = False
                    coords[:] = safe_coords           
              

class AppDMAGMINBH(AppClusterBH):
    overlap_cutoff = None
    genrandom_volume = None
    random_start = False
    
    def __init__(self):       
        AppClusterBH.__init__(self)
        self.quenchParameters={}
        quenchParameters=self.quenchParameters
        quenchParameters["tol"]=1e-4
        quenchParameters["nsteps"]=5000
        quenchParameters["maxErise"]=2e-2
        quenchParameters["maxstep"]=0.1
        quenchParameters["M"]=100
        
        self.quenchRoutine=quenchCrystal
        
    def set_random_start(self, random_start):
        self.random_start = random_start
    
    def initial_coords(self):
        coords = self.create_potential().getCoords()
        if(not self.random_start):
            return coords
        rand = self.create_takestep_reseed()
        rand.takeStep(coords)
        return coords
        
    def create_potential(self):
        return DMACRYSPotential()
    
    def create_takestep_step(self):
        return TakestepDMAGMIN(overlap_cutoff = self.overlap_cutoff, nmols=self.options.displace_nmols)

    def create_takestep_reseed(self):
        return GenRandomCrystal(CoordsAdapter(nrigid=GMIN.getNRigidBody(), nlattice=6, coords=None),
                                volume=self.genrandom_volume,
                                overlap = self.overlap_cutoff)
        
    def create_takestep(self):
        step = self.create_takestep_step()
        if(not self.options.reseed is None):
            reseed = self.create_takestep_reseed()
            return takestep.Reseeding(step, reseed, self.options.reseed)
        return step
    
    def create_basinhopping(self):
        GMIN.initialize()    
        opt = AppClusterBH.create_basinhopping(self)
        self.database.compareMinima = compareMinima
        opt.insert_rejected = True
        addColdFusionCheck(opt)
        return opt
    
    def add_options(self):
        AppBasinHopping.add_options(self)
        self.add_option("-r","--rotate",type="float",
                           dest="rotate", default=3.1415,
                           help="size of rotational displacement",
                           group="Takestep")
        self.add_option("-d","--displace",type="float",
                           dest="displace", default=3.1415,
                           help="carthesian displacement for molecules",
                           group="Takestep")
        self.add_option("--nmols",type="float", dest="displace_nmols", group="Takestep", 
                           help="numper of molecules to displace in each step, default is all")
        self.add_option("--reseed", type="int",dest="reseed",
                           help="give number of maxnoimporve steps to enable reseeding",
                          group="Takestep")
