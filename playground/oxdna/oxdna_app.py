import time
import numpy as np
from pele.application import AppBasinHopping
from pele.potentials import GMINPotential
import oxdnagmin_ as GMIN
from pele import takestep
from math import pi
from pele.utils.rbtools import CoordsAdapter
from pele.utils import rotations
from pele.systems.oxdna import *
from pele.optimize import mylbfgs 

EDIFF=0.01
t0=time.clock()

# this is the base application which controls the program flow
class AppOXDNA(AppBasinHopping):
    #target=-19.40126
    target=None
    default_temperature = 0.4
    
    def __init__(self, *args, **kwargs):
        AppBasinHopping.__init__(self, *args, **kwargs)
        self.quenchRoutine = mylbfgs
        self.potential = GMINPotential(GMIN)
        
        self.quenchParameters["tol"]=1e-4
        self.quenchParameters["M"]=80
        self.quenchParameters["maxErise"]=0.1
        self.quenchRoutine=mylbfgs

        
    # create potentia which calls GMIN
    def create_potential(self):
        return self.potential
    
    # This routine implements how basin hopping should do the step taking
    # Block moves should go in here as well
    def create_takestep(self):    
        # we combine a normal step taking
        # step = OXDNATakestep(displace=self.options.displace, rotate=self.options.rotate)
        group = takestep.BlockMoves()

        #step1 = takestep.AdaptiveStepsize(OXDNATakestep(displace=self.options.displace, rotate=0.), frequency=50)
        step2 = takestep.AdaptiveStepsize(OXDNATakestep(displace=0., rotate=self.options.rotate), frequency=50)
        
        step1 = takestep.AdaptiveStepsize(OXDNAScrewStep(rotate_backbone=0.5, rotate_base=3.))
        group.addBlock(self.options.block1, step1)
        group.addBlock(100, step2)

        # with a generate random configuration
        genrandom = OXDNAReseed()
        # in a reseeding takestep procedure
        reseed = takestep.Reseeding(step2, genrandom, maxnoimprove=self.options.reseed)
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
                  dest="displace", default=.5,
                  help="dispacement for random step taking",
                  group="OXDNA specific")
        self.add_option("-r","--rotate",type="float",
                  dest="rotate", default=0.5,
                  help="rotations for random step taking",
                  group="OXDNA specific")
        self.add_option("--reseed",type="int",
                  dest="reseed", default=1000,
                  help="maximum number of steps with no improvement",
                  group="Basinhopping")
        self.add_option("--block1",type="int",
                  dest="block1", default=100,
                  help="maximum number of steps with no improvement",
                  group="Basinhopping")
        
    def create_basinhopping(self, add_minimum=None, potential=None):
        if(self.target != None):
            add_minimum = self.check_converged

        self.opt = AppBasinHopping.create_basinhopping(self, add_minimum=add_minimum)
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
