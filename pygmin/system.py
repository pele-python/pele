from optparse import OptionParser
from pygmin.storage.database import Database
from pygmin.basinhopping import BasinHopping
from pygmin import takestep

class System(object):
    class Options: 
        bh_temperature = 1.0
        bh_accuracy = 1.e-2
        database="storage.sqlite"
    
    opts = Options()
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
        
        opts = self.opts
        
        print "Initial energy ", pot.getEnergy(coords)

        
        self.database = Database(db=opts.database, accuracy=opts.bh_accuracy, compareMinima=self.compareStructures)

        return BasinHopping(coords, pot, takeStep=step, 
                           temperature=opts.temperature, storage=self.database.minimum_adder())
    
    def parse_options(self, options, args):
        self.opts.temperature = options.temperature


def run_basinhopping(system):
    parser = OptionParser()
    parser.add_option("-T","--temperature",type="float",
                      dest="temperature", default=1.0,
                      help="temperature for the basin hopping run")
    parser.add_option("-n","--nsteps",type="int",
                      dest="nsteps", default=100,
                      help="numper of steps for basin hopping run")
        
    (options, args) = parser.parse_args()
    
    system.parse_options(options, args)
    
    opt = system.create_basinhopping()
    opt.run(options.nsteps)
    