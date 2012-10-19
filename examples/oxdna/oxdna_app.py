import time
import numpy as np
from pygmin.application import AppBasinHopping
from pygmin import defaults
from pygmin.potentials import GMINPotential
import oxdnagmin_ as GMIN
from pygmin import takestep
from math import pi
from pygmin.utils.rbtools import CoordsAdapter
from pygmin.utils import rotations

EDIFF=0.01
t0=time.clock()

# This is the takestep routine for OXDNA. It is a standard rigid body takestep
# routine, but I put it here to be able to start modifying it
class OXDNATakestep(takestep.TakestepInterface):
    def __init__(self, displace=1.0, rotate=0.5*pi):
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

# this is the base application which controls the program flow
class AppOXDNA(AppBasinHopping):
    target=-17.7371736249

    def __init__(self, *args, **kwargs):
        AppBasinHopping.__init__(self, *args, **kwargs)
        self.quenchRoutine = defaults.quenchRoutine
        self.potential = GMINPotential(GMIN)
        
    # create potentia which calls GMIN
    def create_potential(self):
        return self.potential
    
    # This routine implements how basin hopping should do the step taking
    # Block moves should go in here as well
    def create_takestep(self):    
        # we combine a normal step taking
        # step = OXDNATakestep(displace=self.options.displace, rotate=self.options.rotate)
        group = takestep.BlockMoves()

        step1 = takestep.AdaptiveStepsize(OXDNATakestep(displace=self.options.displace, rotate=0.), frequency=50)
        step2 = takestep.AdaptiveStepsize(OXDNATakestep(displace=0., rotate=self.options.rotate), frequency=50)
        group.addBlock(100, step1)
        group.addBlock(100, step2)

        # with a generate random configuration
        genrandom = OXDNAReseed()
        # in a reseeding takestep procedure
        reseed = takestep.Reseeding(group, genrandom, maxnoimprove=self.options.reseed)
        return reseed
    
    # generate some initial coordinates for basin hopping
    def initial_coords(self):
        # get the initial coordinates directly from GMIN
        coords = self.create_potential().getCoords()
        # genearte a random configuration
        OXDNAReseed().takeStep(coords)
        return coords

    # we define some options here which can be set from command line
    # the results will be stored in 
    # self.options.<optionname>
    def add_options(self):
        AppBasinHopping.add_options(self)
        self.add_option("-d","--displace",type="float",
                  dest="displace", default=1.0,
                  help="dispacement for random step taking",
                  group="OXDNA specific")
        self.add_option("-r","--rotate",type="float",
                  dest="rotate", default=pi,
                  help="rotations for random step taking",
                  group="OXDNA specific")
        self.add_option("--reseed",type="int",
                  dest="reseed", default=1000,
                  help="maximum number of steps with no improvement",
                  group="Basinhopping")
        
    def create_basinhopping(self, add_minimum=None):
        if(self.target != None):
            add_minimum = self.check_converged

        self.opt = AppBasinHopping.create_basinhopping(self, add_minimum=self.check_converged)
	return self.opt
 
    def check_converged(self, E, coords):
        if(E<(self.target+EDIFF)):
            print "#found minimum"
            t1= time.clock()
            timespent= t1 - t0
            print "quenches functioncalls timespent"
            print "%d %d %f"%(self.opt.stepnum, self.potential.ncalls, timespent)
            exit()


        
if __name__ == "__main__":
    # let GMIN do some initialization
    GMIN.initialize()

    # create a OXDNA application
    app = AppOXDNA()
    # and run it
    app.execute()