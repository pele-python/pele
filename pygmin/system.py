from optparse import OptionParser
from pygmin.storage.database import Database
from pygmin.basinhopping import BasinHopping
from pygmin import takestep

class System(object):   
    options = None
    accuracy = 1e-3
    compareStructures = None
    database = None
    
    
    def create_potential(self):
        raise Exception("system does not implement create_potential")
        
    def create_takestep(self):
        return takestep.RandomDisplacement()
    
    def initial_coords(self):
        raise Exception("system does not implement initial_coords")

    def create_basinhopping(self):
        pot = self.create_potential()
        coords = self.initial_coords()
        step = self.create_takestep()
        
        opts = self.options
        
        print "Initial energy ", pot.getEnergy(coords)

        
        self.database = Database(db=opts.database, accuracy=self.accuracy, compareMinima=self.compareStructures)

        return BasinHopping(coords, pot, takeStep=step, 
                           temperature=opts.temperature, storage=self.database.minimum_adder())
    
    def parse_options(self, options, args):
        self.options = options
        
class ClusterSystem(System):
    def create_takestep(self):
        step = takestep.UniformDisplacement(stepsize=self.options.stepsize)
        if(not self.options.reseed is None):
            reseed = takestep.RandomCluster(volume=1.0)
            return takestep.Reseeding(step, reseed, self.options.reseed)
        return step
        
    def add_basinhopping_options(self, parser, bh_opts):
        ts_opts = parser.add_option_group("Takestep")        
        ts_opts.add_option("-s","--step-size",type="float",
                           dest="stepsize", default=1.0,
                           help="takestep stepsize")
        ts_opts.add_option("--reseed", type="int",dest="reseed",
                           help="give number of maxnoimporve steps to enable reseeding")

def run_basinhopping(system):
    parser = OptionParser()
    parser.add_option("--db",type="string",
                      dest="database", default="storage.sqlite",
                      help="database to store results")

    bh_opts = parser.add_option_group("Basin Hopping")
    bh_opts.add_option("-T","--temperature",type="float",
                      dest="temperature", default=1.0,
                      help="temperature for the basin hopping run")
    bh_opts.add_option("-n","--nsteps",type="int",
                      dest="nsteps", default=100,
                      help="numper of steps for basin hopping run")
    
    system.add_basinhopping_options(parser, bh_opts)
        
    (options, args) = parser.parse_args()
    
    system.parse_options(options, args)
    
    opt = system.create_basinhopping()
    opt.run(options.nsteps)
    