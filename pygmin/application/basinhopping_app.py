from pygmin.basinhopping import BasinHopping
from pygmin import takestep
from application_base import Application
from pygmin.storage.database import Database
from pygmin import defaults

class AppBasinHopping(Application):
    accuracy = 1e-3
    quenchParameters=defaults.quenchParams
    quenchRoutine=defaults.quenchRoutine
    
    def create_takestep(self):
        return takestep.RandomDisplacement()
    
    def initial_coords(self):
        raise Exception("system does not implement initial_coords")

    def create_basinhopping(self, add_minimum=None):
        pot = self.create_potential()
        coords = self.initial_coords()
        step = self.create_takestep()
        
        opts = self.options
        
        print "Initial energy ", pot.getEnergy(coords)

        if(add_minimum is None):
            add_minimum = self.database.minimum_adder()
        return BasinHopping(coords, pot, takeStep=step, 
                           temperature=opts.temperature, storage=add_minimum,
                           quenchRoutine=self.quenchRoutine,
                           quenchParameters = self.quenchParameters)
        
    def add_options(self):
        self.add_option("--db",type="string",
                          dest="database", default="storage.sqlite",
                          help="database to store results")
    
        self.add_option("-T","--temperature",type="float",
                          dest="temperature", default=1.0,
                          help="temperature for the basin hopping run",
                          group="Basin Hopping")
        self.add_option("-n","--nsteps",type="int",
                          dest="nsteps", default=100,
                          help="numper of steps for basin hopping run")
                
    def run(self):        
        opt = self.create_basinhopping()        
        opt.run(self.options.nsteps)
        
        
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
