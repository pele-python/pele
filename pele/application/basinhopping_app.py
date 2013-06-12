from pele.basinhopping import BasinHopping
from pele import takestep
from application_base import Application
from pele.storage.database import Database
from pele import defaults
import time

class AppBasinHopping(Application):
    accuracy = 1e-3
    default_temperature = 1.0
    target_energy = None
    
    def __init__(self):
        self.quenchParameters=defaults.quenchParams
        self.quenchRoutine=defaults.quenchRoutine
        
    def create_takestep(self):
        return takestep.RandomDisplacement()
    
    def initial_coords(self):
        raise Exception("system does not implement initial_coords")

    def create_basinhopping(self, add_minimum=None, potential=None):
        if(potential is None):
            potential = self.create_potential()
            
        coords = self.initial_coords()
        step = self.create_takestep()
        
        opts = self.options
        
        print "Initial energy ", potential.getEnergy(coords)
        print "Creating basin hopping with quencher"
        print "------------------------------------"
        print self.quenchRoutine
        for key, value in self.quenchParameters.iteritems():
            print key, "=", value         
        print "------------------------------------"
        
        if self.target_energy is not None:
            add_minimum = self.check_converged
        
        if(add_minimum is None):
            add_minimum = self.database.minimum_adder()
        
        return BasinHopping(coords, potential, takeStep=step, 
                           temperature=opts.temperature, storage=add_minimum,
                           quenchRoutine=self.quenchRoutine,
                           quenchParameters = self.quenchParameters)
        
    def add_options(self):
        self.add_option("--db",type="string",
                          dest="database", default="storage.sqlite",
                          help="database to store results")
    
        self.add_option("-T","--temperature",type="float",
                          dest="temperature", default=self.default_temperature,
                          help="temperature for the basin hopping run",
                          group="Basin Hopping")
        self.add_option("-n","--nsteps",type="int",
                          dest="nsteps", default=1000,
                          help="numper of steps for basin hopping run")
        
    def check_converged(self, E, coords):
        if(E<(self.target_energy)):
            print "#found minimum"
            t1= time.time()
            timespent= t1 - self.t0
            print "quenches functioncalls timespent"
            print "%d %d %f"%(self.opt.stepnum, self.potential.ncalls, timespent)
            exit()
                
    def run(self):        
        self.t0 = time.time()
        self.potential = self.create_potential()
        self.opt = self.create_basinhopping(potential=self.potential)        
        self.opt.run(self.options.nsteps)
        
        
class AppClusterBH(AppBasinHopping):
    def create_takestep_step(self):
        return takestep.UniformDisplacement(stepsize=self.options.stepsize) 
    def create_takestep_reseed(self):
        return takestep.RandomCluster(volume=1.0)
        
    def create_takestep(self):
        step = self.create_takestep_step()
        if(not self.options.reseed is None):
            reseed = self.create_takestep_reseed()
            return takestep.Reseeding(step, reseed, self.options.reseed)
        return step
        
    def add_options(self):
        AppBasinHopping.add_options(self)
        self.add_option("-s","--step-size",type="float",
                           dest="stepsize", default=1.0,
                           help="takestep stepsize",
                           group="Takestep")
        self.add_option("--reseed", type="int",dest="reseed",
                           help="give number of maxnoimporve steps to enable reseeding",
                          group="Takestep")
